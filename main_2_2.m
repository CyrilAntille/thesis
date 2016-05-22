%% 2.2: Motion between beams
% addpath '/uio/hume/student-u04/cyrila/Documents/MATLAB/MasterThesis/Field_II_ver_3_24' -end
addpath 'C:\Users\Cyril\Documents\MATLAB\Field_II_ver_3_24' -end

save_all_data = false; % Can take multiple GB of memory
enable_plots = true;
grayscale_plots = true;

bfm = [0,1,3,4]; % 0=DAS, 1=MV, 2=MV-Multibeam, 3=IAA, 4=IAA_500pts
methods_set = {'DAS','MV','IAA','IAA 500pts'}; % Must correspond to bfm

disable_multiprocess = false; % Must disable automatic creation of 
% parallel pool for parfor (in parallel preferences).

%% 0. Load data
% If raw/DA/BF data exists, the script won't recreate the data.
% -> Can just load data from file if available

% load 2_1_speckle.mat % Loads raw speckle data
% bf_data_files = {'2_1_61_91_noSpeckle.mat', '2_1_101_171_noSpeckle.mat', ...
%     '2_1_181_211_noSpeckle.mat'};
% [ shift, num_beams, data_phantoms, data_DA, data_BF ] ...
%     = mergeData(bf_data_files, false); % Loads and merges DA(S) (and BF) data

%% 1. Create Speckle raw data
field_init(0);

create_speckle = false; 
if create_speckle && exist('speckle_raw', 'var') ~= 1
    save_speckle = true;
    fprintf('Creating raw speckle image. This can take many hours if many points.\n')
    P = Parameters();
    P.Seed = 2;
    P.NumPoints = 10^6;
    speckle_phantom = PlanewaveVesselPhantom(P, 0, P.NumPoints, P.Seed);
    speckle_raw = CalcRespAll(P, speckle_phantom);
    
    if save_speckle
        output_file = ['2_1_speckle_' num2str(P.Seed) '_10-' ...
            num2str(log10(P.NumPoints)) '.mat'];
        fprintf('\nNSaving speckle raw data to file.')
        save(output_file, 'P', 'speckle_raw', 'speckle_phantom', '-v7.3')
    end
end
fprintf('\n============================================================\n')

add_speckle = true;
if exist('speckle_raw', 'var') ~= 1 || isempty(speckle_raw)
    P = Parameters();
    add_speckle = false;
    speckle_raw = [];
end

fprintf('Initializing Field II for all workers.\n')
if disable_multiprocess
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
if exist('data_DA', 'var') ~= 1
    fprintf('Creating raw and DA(S) images. This can take hours if many beam setups (NThetas).\n')
    shift = Shift(ShiftType.LateralCst, 0, -1);
    if exist('num_beams', 'var') ~= 1
        num_beams = 61:10:71;
    end
    
    pts_theta = [0, 0]; % Add a theta (in degrees) for each point
    pts_range = [P.Tx.FocRad, 50*1e-3]; % Add a range (in m) for each point
    scat_pts = zeros([length(pts_theta) 3]);
    for pidx = 1:length(pts_theta)
        pt = [sind(pts_theta(pidx)) 0 cosd(pts_theta(pidx))] ...
            * pts_range(pidx);
        scat_pts(pidx, :) = pt;
    end
    original_phantom = PointPhantom(scat_pts, 1+db2mag(20));
    % -> pts 20dB over speckle
    
    data_phantoms = original_phantom;
    data_DA = cell([1, length(num_beams)]);
    parfor b=1:length(num_beams)
        Pb = copyStruct(P);
        Pb.Tx.NTheta = num_beams(b);
        Pb.Tx.SinTheta = linspace(-Pb.Tx.SinThMax, Pb.Tx.SinThMax, Pb.Tx.NTheta);
        Pb.Tx.Theta = asin(Pb.Tx.SinTheta);
        fprintf('\nNTheta: %d.\n', Pb.Tx.NTheta)
        nstart = tic;
        
        b_shift = Shift(shift.type, shift.val, Pb.Tx.NTheta);
        b_DA = cell([1, b_shift.num_shifts]);
        shifts = b_shift.getShifts(Pb);
        for s=1:b_shift.num_shifts
            Ps = copyStruct(Pb);
            Ps.Tx.NTheta = 1;
            Ps.Tx.SinTheta = - Pb.Tx.SinThMax + shifts(s);
            Ps.Tx.Theta = asin(Ps.Tx.SinTheta);
            s_phantom = shift.shiftPositions(original_phantom, shifts(s));
            s_raw = CalcRespAll(Ps, s_phantom);
            if add_speckle
                img = speckle_raw.image;
                % Assumes scatterers in speckle -> s_raw.image smaller than img.
                img(1:size(s_raw.image, 1),:) = img(1:size(s_raw.image, 1),:) + s_raw.image;
                s_raw.image = img;
            end
            b_DA{s} = BeamformAll(Ps, s_raw);
        end
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
        data_DA{b} = n_DA;
        nend = toc(nstart);
        fprintf('\nNTheta: %d. Time: %d minutes and %f seconds', ...
            Pb.Tx.NTheta, floor(nend/60), rem(nend,60))
    end
    if save_all_data
        speckle_name = '_noSpeckle';
        if add_speckle
            speckle_name = ['_speckle' num2str(P.Seed)];
        end
        output_file = ['2_1_' num2str(num_beams(1)) '_' ...
            num2str(num_beams(end)) speckle_name '.mat'];
        fprintf('\nNSaving DA(S) data to file.')
        save(output_file, 'shift', 'num_beams', 'data_phantoms', ...
            'data_DA', '-v7.3')
    end
end
field_end();
parfor w=1:pool.NumWorkers
    field_end();
end
fprintf('\n============================================================\n')

%% 3. Beamform data: DAS, MV, IAA
if exist('data_BF', 'var') ~= 1
    fprintf('Beamforming all DA(S) images.\n')
    data_BF = cell([1, length(bfm)]);
    for m=1:length(bfm)
        bf_method = bfm(m);
        fprintf('----------------------\n')
        fprintf('Beamforming method: %s.\n', methods_set{m})
        mstart = tic;
        
        m_BF = cell([1, length(num_beams)]);
        parfor b=1:length(num_beams)
            Pb = copyStruct(P);
            Pb.Tx.NTheta = num_beams(b);
            Pb.Tx.SinTheta = linspace(-Pb.Tx.SinThMax,Pb.Tx.SinThMax,Pb.Tx.NTheta);
            Pb.Tx.Theta = asin(Pb.Tx.SinTheta);
            b_DA = data_DA{b};
            m_BF{b} = ComputeBF(b_DA.image, Pb, bf_method);
        end
        data_BF{m} = m_BF;
        mend = toc(mstart);
        fprintf('%s: %d minutes and %f seconds\n', methods_set{m}, ...
            floor(mend/60), rem(mend,60))
    end
    if save_all_data
        speckle_name = '_noSpeckle';
        if add_speckle
            speckle_name = ['_speckle' num2str(P.Seed)];
        end
        output_file = ['2_1_' num2str(num_beams(1)) '_' ...
            num2str(num_beams(end)) speckle_name '.mat'];
        fprintf('\nNSaving beamformed data to file.')
        save(output_file, 'data_BF', '-append')
    end
end
fprintf('\n============================================================\n')

%% Plots
if enable_plots
    % BF images
    figure;
    for m=1:length(methods_set)
        m_BF = data_BF{m};
        for b=1:length(num_beams)
            Pb = copyStruct(P);
            Pb.Tx.NTheta = num_beams(b);
            Pb.Tx.SinTheta = linspace(-Pb.Tx.SinThMax,Pb.Tx.SinThMax,Pb.Tx.NTheta);
            Pb.Tx.Theta = asin(Pb.Tx.SinTheta);

            b_BF = m_BF{b};
            b_DA = data_DA{b};

            thetaRange = Pb.Tx.Theta;
            if strcmp(methods_set(m),'IAA 500pts')
                thetaRange = linspace(Pb.Tx.Theta(1), Pb.Tx.Theta(end), 500);
            end
            warning('off')
            [scanConvertedImage, Xs, Zs] = getScanConvertedImage(b_BF, ...
                thetaRange, 1e3 * b_DA.Radius, 2024, 2024, 'spline');
            warning('on')
            img = db(abs(scanConvertedImage));
            imagesc(Xs, Zs, img)
            xlabel('azimuth [mm]');

%                 imagesc(rad2deg(thetaRange),1e3 * s_DA.Radius, db(abs(s_BF)));
%                 xlim([-15 15])
%                 ylim([35 45])
%                 xlabel('angle [deg]');

            caxis([-130  -70]);
            colorbar
            if grayscale_plots
                colormap(gray)
            else 
                colormap(jet)
            end 
            ylabel('range [mm]');
            title(['Method :', methods_set{m}, ', No Beams: ', ...
                int2str(num_beams(b))])
            pause
        end
    end
end