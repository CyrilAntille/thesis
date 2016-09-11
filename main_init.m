%% 3.0: Parameters initialization for all experiments
addpath ../Field_II_ver_3_24/ -end

if ~exist('mainP', 'var')
    mainP = MainParameters();
end

%% 1. FieldII initialization and Speckle raw data creation
field_init(0);

if mainP.speckle_create && ~exist('speckle_raw', 'var')
    fprintf('Creating raw speckle image. \n')
    fprintf('Speckle can take a long time to simulate.\n')
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
    load(mainP.speckle_file)
end
fprintf('\n============================================================\n')

if ~exist('speckle_raw', 'var')
    speckle_raw_image = [];
else
    speckle_raw_image = speckle_raw.image;
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
parfor w=1:pool.NumWorkers
    field_init(0);
end
fprintf('\n============================================================\n')

%% 2. Create raw and DA(S) data with point scatterers
if ~exist('data_DA', 'var')
    fprintf('Creating raw and DA(S) images. This can take hours if many beam setups (NThetas).\n')
    scat_pts = zeros([length(mainP.pts_range) 3]);
    pts_angle = 0; % Must be 0 to guaranty a transmit beam hitting the pt
    for pidx = 1:length(mainP.pts_range)
        scat_pts(pidx, :) = [sind(pts_angle) 0 cosd(pts_angle)] ...
            * mainP.pts_range(pidx) * 1e-3;
    end
    original_phantom = PointPhantom(scat_pts, 1+db2mag(30));
    % -> pts 30dB over speckle
    
    if mainP.shift_per_beam
        data_phantoms = original_phantom;
    else
        data_phantoms = cell([1, length(mainP.num_beams)]);
    end
    data_DA = cell([1, length(mainP.num_beams)]);
    for b=1:length(mainP.num_beams)
        Pb = mainP.copyP(mainP.num_beams(b));
        fprintf('\nNTheta: %d.\n', Pb.Tx.NTheta)
        nstart = tic;
        
        if mainP.shift_per_beam
            b_shift = Shift(mainP.shift.type, mainP.shift.val, ...
                Pb.Tx.NTheta, mainP.shift.direction);
            beam_shift = Pb.Tx.SinTheta(2) - Pb.Tx.SinTheta(1);
        else
            b_shift = mainP.shift;
            b_phantoms = cell([1, b_shift.num_shifts]);
        end
        
        b_DA = cell([1, b_shift.num_shifts]);
        shifts = b_shift.getShifts(Pb);
        parfor s=1:b_shift.num_shifts
            % Raw
            s_phantom = b_shift.shiftPositions(original_phantom, shifts(s));
            s_raw = CalcRespAll(Pb, s_phantom);
            if ~isempty(speckle_raw_image)
                new_image = speckle_raw_image;
                new_image(1:size(s_raw.image),:) = ...
                    new_image(1:size(s_raw.image),:) + s_raw.image;
                s_raw.image = new_image;
                %TODO: Maybe always include focus range to image!
            end
            if mainP.shift_per_beam
                Ps = mainP.copyP(1);
                Ps.Tx.SinTheta = - Ps.Tx.SinThMax + beam_shift * (s-1);
                Ps.Tx.Theta = asin(Ps.Tx.SinTheta);
            else
                Ps = Pb;
                b_phantoms{s} = s_phantom;
            end
            % DA(S)
            b_DA{s} = BeamformAll(Ps, s_raw);
        end
        if mainP.shift_per_beam
            % Need to merge beamformed images into a single one.
            n_DA = b_DA{1};
            for s=2:b_shift.num_shifts
                radius_length = length(b_DA{s}.Radius); % This number may vary
                if length(n_DA.Radius) > radius_length
                    n_DA.image = n_DA.image(:,1:radius_length,:);
                    n_DA.Radius = b_DA{s}.Radius;
                elseif radius_length > length(n_DA.Radius)
                    b_DA{s}.image = b_DA{s}.image(:,1:length(n_DA.Radius),:);
                end
                n_DA.image = vertcat(n_DA.image, b_DA{s}.image);
            end
            data_DA{b} = n_DA; % single image per num_beam value
        else
            data_DA{b} = b_DA; % num_shifts images per num_beam value
            data_phantoms{b} = b_phantoms;
        end
        nend = toc(nstart);
        fprintf('\nNTheta: %d. Time: %d minutes and %f seconds', ...
            Pb.Tx.NTheta, floor(nend/60), rem(nend,60))
    end
    if mainP.save_all_data
        output_file = mainP.outputFileName(~isempty(speckle_raw_image));
        fprintf('\nNSaving DA(S) data  into: %s\n', output_file)
        save(output_file, 'mainP', 'shift', 'num_beams', 'data_phantoms',...
            'data_DA', '-v7.3')
    end
end
field_end();
parfor w=1:pool.NumWorkers
    field_end();
end
fprintf('\n============================================================\n')

%% 3. Beamform data: DAS, MV, IAA
if ~exist('data_BF', 'var')
    fprintf('Beamforming all DA(S) images.\n')
    data_BF = cell([1, length(mainP.methods_set)]);
    for m=1:length(mainP.methods_set)
        bf_method = mainP.methods_set{m};
        fprintf('----------------------\n')
        fprintf('Beamforming method: %s.\n', bf_method)
        mstart = tic;
        
        m_BF = cell([1, length(mainP.num_beams)]);
        parfor b=1:length(mainP.num_beams)
            mainPs = copyStruct(mainP);
            mainPs.P = mainPs.copyP(mainPs.num_beams(b));
            b_DA = data_DA{b};
            if mainP.shift_per_beam
                b_BF = ComputeBF(b_DA.image, mainPs, bf_method);
            else
                b_BF = cell([1, mainPs.shift.num_shifts]);
                for s=1:mainPs.shift.num_shifts
                    b_BF{s} = ComputeBF(b_DA{s}.image, mainPs, bf_method);
                end
            end
            m_BF{b} = b_BF;
        end
        data_BF{m} = m_BF;
        mend = toc(mstart);
        fprintf('%s: %d minutes and %f seconds\n', mainP.methods_set{m},...
            floor(mend/60), rem(mend,60))
    end
    if mainP.save_all_data
        output_file = mainP.outputFileName(~isempty(speckle_raw_image));
        fprintf('\nNSaving beamformed data into: %s\n', output_file)
        save(output_file, 'data_BF', '-append')
    end
end
fprintf('\n============================================================\n')
fprintf('main_init finished!\n')

