%% 3.0: Parameters initialization for all experiments
addpath ../Field_II_ver_3_24/ -end

if ~exist('mainP', 'var')
    mainP = MainParameters();
end

%% 1. FieldII initialization and Speckle raw data creation
fprintf('============================================================\n')
fprintf('Initializing Field II for main process\n')
evalc('field_init(0)');

if mainP.speckle_create && ~exist('speckle_raw', 'var')
    fprintf('Creating raw speckle image\n')
    fprintf('Speckle can take a long time to simulate\n')
    speckle_phantom = PlanewaveVesselPhantom(mainP.P, 0, ...
        mainP.speckle_numPts, mainP.speckle_seed);
    speckle_raw = CalcRespAll(mainP.P, speckle_phantom);
    
    if mainP.speckle_save
        output_file = ['2_1_speckle_' num2str(mainP.speckle_seed) '_10-'...
            num2str(log10(mainP.speckle_numPts)) '.mat'];
        output_file = strcat(mainP.save_folder, output_file);
        fprintf('\nNSaving speckle raw data into: %s\n', output_file)
        save(output_file, 'mainP', 'speckle_raw', 'speckle_phantom', '-v7.3')
    end
elseif mainP.speckle_load && ~exist('speckle_raw', 'var')
    load(mainP.speckle_file, 'speckle_raw')
end

if exist('speckle_raw', 'var')
    resh_speckle = reshapeRawImg(mainP, speckle_raw);
    speckle_raw_image = resh_speckle.image;
else
    mock_speckle = struct; mock_speckle.times = 0;
    mock_speckle.image = []; mock_speckle.LConvs = 0;
    mock_speckle = reshapeRawImg(mainP, mock_speckle);
    speckle_raw_image = mock_speckle.image; % zero valued image
end
% speckle_raw_image is created to avoid copying the whole speckle_raw
% structure when running parallel loops.

fprintf('Initializing Field II for all workers.\n')
if mainP.disable_multiprocess
    pool = struct; pool.NumWorkers = 1;
else
    if isempty(gcp('nocreate'))
        parpool;
    end
    pool = gcp;
end
f = parfevalOnAll(pool,@field_init,0,0);
wait(f);
fprintf('============================================================\n')

%% 2. Create raw and DA(S) data with point scatterers
if ~exist('data_DA', 'var')
    fprintf('Creating raw and DA(S) images\n')
    nstart = tic;
    % Scatterer points creation
    scat_pts = zeros([length(mainP.pts_range) 3]);
    for pidx = 1:length(mainP.pts_range)
        scat_pts(pidx, :) = [mainP.pts_azimuth(pidx) 0 ...
            mainP.pts_range(pidx)] .* 1e-3;
    end
    phantom = PointPhantom(scat_pts, 1+db2mag(mainP.pts_gain));
    % -> pts mainP.pts_gain over speckle

    fprintf('Progress:\n');
    fprintf([repmat('.', 1, mainP.shift.num_shifts) '\n\n']);
    shifts = mainP.shift.getShifts(mainP.P);
    data_DA = cell([1, mainP.shift.num_shifts]);
    if mainP.shift_per_beam
        data_phantom = phantom;
    else
        data_phantom = cell([1, mainP.shift.num_shifts]);
    end
    parfor s=1:mainP.shift.num_shifts
        % Raw
        s_phantom = mainP.shift.shiftPositions(phantom, shifts(s));
        if ~mainP.shift_per_beam
            data_phantom{s} = s_phantom;
        end
        s_raw = CalcRespAll(mainP.P, s_phantom);
        s_raw = reshapeRawImg(mainP, s_raw);
        s_raw.image = s_raw.image + speckle_raw_image;
        if mainP.shift_per_beam
            beam_shift = mainP.P.Tx.SinTheta(2) - mainP.P.Tx.SinTheta(1);
            Ps = copyStruct(mainP.P); Ps.Tx.NTheta = 1;
            Ps.Tx.SinTheta = - Ps.Tx.SinThMax + beam_shift * (s-1);
            Ps.Tx.Theta = asin(Ps.Tx.SinTheta);
        else
            Ps = mainP.P;
        end
        % DA(S)
        data_DA{s} = BeamformAll(Ps, s_raw);
        fprintf('\b|\n');
    end
    if mainP.shift_per_beam
        % Need to merge beamformed images into a single one.
        n_DA = data_DA{1};
        for s=2:mainP.shift.num_shifts
            n_DA.image = vertcat(n_DA.image, data_DA{s}.image);
        end
        data_DA = n_DA; % single image
    end
    nend = toc(nstart);
    fprintf('Raw and DA(s) creation: Time: %d minutes and %0.3f seconds\n', ...
        floor(nend/60), rem(nend,60))
    if mainP.save_all_data
        output_file = mainP.outputFileName(false);
        fprintf('\nNSaving DA(S) data  into: %s\n', output_file)
        save(output_file, 'mainP', 'data_phantom', 'data_DA', '-v7.3')
    end
end
field_end();
f = parfevalOnAll(pool,@field_end,0);
wait(f);
fprintf('============================================================\n')

%% 3. Beamform data: DAS, MV, IAA
if ~exist('data_BF', 'var')
    fprintf('Beamforming all DA(S) images\n')
    data_BF = cell([1, length(mainP.methods_set)]);
    parfor m=1:length(mainP.methods_set)
        bf_method = mainP.methods_set{m};
%         fprintf('----------------------\n')
%         fprintf('Beamforming method: %s.\n', bf_method)
        mstart = tic;
        if mainP.shift_per_beam
            m_BF = ComputeBF(data_DA.image, mainP, bf_method, 0);
        else
            m_BF = cell([1, mainP.shift.num_shifts]);
            for s=1:mainP.shift.num_shifts
                m_BF{s} = ComputeBF(data_DA{s}.image, mainP, bf_method, 0);
            end
        end
        data_BF{m} = m_BF;
        mend = toc(mstart);
        fprintf('%s: %d minutes and %0.3f seconds\n', mainP.methods_set{m},...
            floor(mend/60), rem(mend,60))
    end
end
%% Test
% prefix = mainP.files_prefix;
% mainP.files_prefix = strcat('before_', mainP.files_prefix);
% upval = mainP.interp_upsample;
% mainP.interp_upsample = 0;
% plotBFImages(mainP, data_DA, data_BF)
% mainP.interp_upsample = upval;
% mainP.files_prefix = prefix;
%%
if true
    fprintf('Images upsampling and normalization\n')
    % mainP.norm_variant: 1 = Standalone normalization, 2 = DAS-Based normalization,
    % 3 = First-shift based normalization (only if shift_per_beam = false)
    if mainP.shift_per_beam
        norm_info = [max(data_BF{1}(:)), 0, mainP.normalize_bfim];
        % norm_info = [orig_max, norm_factor, do_normalization]
        data_peaks = cell([1, mainP.num_beams]);
        for m=1:length(mainP.methods_set)
            if mainP.norm_variant == 2 && m > 1
                mainP.normalize_bfim = false;
                m_BF = normalizeBFImage(mainP, data_BF{m}, ...
                    data_DA.Radius, data_phantom);
                if norm_info(3)
                    m_BF = m_BF .* norm_info(2);
%                     m_BF = m_BF - norm_info(2); % dB
                end
                mainP.normalize_bfim = norm_info(3);
            else
                fprintf('whaaat - ')
                m_BF = normalizeBFImage(mainP, data_BF{m}, ...
                    data_DA.Radius, data_phantom);
                norm_info(2) = max(m_BF(:)) / norm_info(1);
%                 norm_info(2) = max(m_BF(:)) - norm_info(1); % dB
            end
            data_peaks{m} = computePeaksInfo(mainP, data_phantom, ...
                data_DA.Radius, m_BF);
            data_BF{m} = m_BF;
        end
    else
        norm_info = [max(data_BF{1}{1}(:)), 0, mainP.normalize_bfim];
        % norm_info = [orig_max, norm_factor, do_normalization]
        data_peaks = cell([1, mainP.shift.num_shifts]);
        for m=1:length(mainP.methods_set)
            if mainP.norm_variant == 3
                norm_info(1) = max(data_BF{m}{1}(:));
            end
            for s=1:mainP.shift.num_shifts
                s_BF =  data_BF{m}{s};
                if s > 1 && (mainP.norm_variant == 2 || mainP.norm_variant == 3)
                    mainP.normalize_bfim = false;
                    s_BF = normalizeBFImage(mainP, s_BF, ...
                        data_DA{s}.Radius, data_phantom{s});
                    if norm_info(3)
                        s_BF = s_BF .* norm_info(2);
%                         s_BF = s_BF - norm_info(2); % dB
                    end
                    mainP.normalize_bfim = norm_info(3);
                else
                    s_BF = normalizeBFImage(mainP, s_BF, ...
                        data_DA{s}.Radius, data_phantom{s});
                    norm_info(2) = max(s_BF(:)) / norm_info(1);
%                     norm_info(2) = max(s_BF(:)) - norm_info(1); % dB
                end
                data_peaks{m}{s} = computePeaksInfo(mainP, ...
                    data_phantom{s}, data_DA{s}.Radius, s_BF);
                data_BF{m}{s} = s_BF;
            end
        end
    end
    if mainP.save_all_data
        output_file = mainP.outputFileName(false);
        fprintf('\nNSaving beamformed data into: %s\n', output_file)
        save(output_file, 'data_BF', 'data_peaks', '-append')
    end
end
fprintf('============================================================\n')
fprintf('main_init finished!\n')

