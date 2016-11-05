%% 3.0: Parameters initialization for all experiments
addpath ../Field_II_ver_3_24/ -end

if ~exist('mainP', 'var')
    mainP = MainParameters();
end

%% 1. FieldII initialization and Speckle raw data creation
fprintf('============================================================\n')
fprintf('Initializing Field II for main process\n')
% field_init(0);
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
%     speckle_raw_image = speckle_raw.image;
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
        scat_pts(pidx, :) = [mainP.pts_azimuth(pidx) * 1e-3 0 ...
            mainP.pts_range(pidx) * 1e-3];
    end
    data_phantom = PointPhantom(scat_pts, 1+db2mag(mainP.pts_gain));
    % -> pts mainP.pts_gain over speckle

    fprintf('Progress:\n');
    fprintf([repmat('.', 1, mainP.shift.num_shifts) '\n\n']);
    shifts = mainP.shift.getShifts(mainP.P);
    data_DA = cell([1, mainP.shift.num_shifts]);
    parfor s=1:mainP.shift.num_shifts
        % Raw
        s_phantom = mainP.shift.shiftPositions(data_phantom, shifts(s));
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
            m_BF = normalizeBFImage(m_BF, data_DA.Radius*1000);
        else
            m_BF = cell([1, mainP.shift.num_shifts]);
            for s=1:mainP.shift.num_shifts
                m_BF{s} = ComputeBF(data_DA{s}.image, mainP, bf_method, 0);
                m_BF{s} = normalizeBFImage(m_BF{s}, data_DA{s}.Radius*1000);
            end
        end
        data_BF{m} = m_BF;
        mend = toc(mstart);
        fprintf('%s: %d minutes and %0.3f seconds\n', mainP.methods_set{m},...
            floor(mend/60), rem(mend,60))
    end
    if mainP.save_all_data
        output_file = mainP.outputFileName(false);
        fprintf('\nNSaving beamformed data into: %s\n', output_file)
        save(output_file, 'data_BF', '-append')
    end
end
fprintf('============================================================\n')
fprintf('main_init finished!\n')

