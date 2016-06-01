%% 2.2: Motion between beams
addpath '/uio/hume/student-u04/cyrila/Documents/MATLAB/MasterThesis/Field_II_ver_3_24' -end
% addpath 'C:\Users\Cyril\Documents\MATLAB\Field_II_ver_3_24' -end

if exist('save_all_data', 'var') ~= 1
    save_all_data = false; % Can take multiple GB of memory
end
if exist('enable_plots', 'var') ~= 1
    enable_plots = true;
end

bfm = [0,1,4,5];
% 0=DAS, 1=MV, 2=MV-MB, 3=IAA-MBSB, 3=IAA-MBMB, 4=IAA-MBMB-Upsampled
% methods_set must correspond to bfm.
methods_set = {'DAS','MV', 'IAA-MBMB','IAA-MBMB-Upsampled'};


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

if exist('speckle_raw', 'var') ~= 1
    P = Parameters();
    speckle_raw_image = [];
else
    speckle_raw_image = speckle_raw.image;
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
    if exist('shift', 'var') ~= 1
%         shift = Shift(ShiftType.LateralCst, 5*1e-3 / 101, -1);
        shift = Shift(ShiftType.LateralSpeed, 1, -1);
    end
    if exist('num_beams', 'var') ~= 1
        num_beams = 101;
%         num_beams = [71 111];
    end
    pts_theta = [15, 15]; % Add a theta (in degrees) for each point
    pts_range = [P.Tx.FocRad, 50*1e-3]; % Add a range (in m) for each point
    scat_pts = zeros([length(pts_theta) 3]);
    for pidx = 1:length(pts_theta)
        pt = [sind(pts_theta(pidx)) 0 cosd(pts_theta(pidx))] ...
            * pts_range(pidx);
        scat_pts(pidx, :) = pt;
    end
    original_phantom = PointPhantom(scat_pts, 1+db2mag(30));
    % -> pts 30dB over speckle
    
    data_phantoms = original_phantom;
    data_DA = cell([1, length(num_beams)]);
    for b=1:length(num_beams)
        Pb = copyStruct(P);
        Pb.Tx.NTheta = num_beams(b);
        Pb.Tx.SinTheta = linspace(-Pb.Tx.SinThMax, Pb.Tx.SinThMax, Pb.Tx.NTheta);
        Pb.Tx.Theta = asin(Pb.Tx.SinTheta);
        fprintf('\nNTheta: %d.\n', Pb.Tx.NTheta)
        nstart = tic;
        
        b_shift = Shift(shift.type, shift.val, Pb.Tx.NTheta);
        b_DA = cell([1, b_shift.num_shifts]);
        shifts = b_shift.getShifts(Pb);
        beam_shift = Pb.Tx.SinTheta(2) - Pb.Tx.SinTheta(1);
        parfor s=1:b_shift.num_shifts
            % Raw
            s_phantom = b_shift.shiftPositions(original_phantom, shifts(s));
            s_raw = CalcRespAll(Pb, s_phantom);
            if ~isempty(speckle_raw_image)
                new_image = speckle_raw_image;
                new_image(1:size(s_raw.image),:) = ...
                    new_image(1:size(s_raw.image),:) + s_raw.image;
                s_raw.image = new_image;
            end
            % DA(S)
            Ps = copyStruct(Pb);
            Ps.Tx.NTheta = 1;
            Ps.Tx.SinTheta = - Pb.Tx.SinThMax + beam_shift * (s-1);
            Ps.Tx.Theta = asin(Ps.Tx.SinTheta);
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
    parfor m=1:length(bfm)
        bf_method = bfm(m);
        fprintf('----------------------\n')
        fprintf('Beamforming method: %s.\n', methods_set{m})
        mstart = tic;
        
        m_BF = cell([1, length(num_beams)]);
        for b=1:length(num_beams)
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

%% 3dB width computations and Plots
if enable_plots
    figure('units','normalized','position',[.2 .3 .5 .3])
end
max_peak_DAS = 0;
points_3dBwidth = cell([length(methods_set), length(num_beams)]);
for m=1:length(methods_set)
    m_BF = data_BF{m};
    for b=1:length(num_beams)
        Pb = copyStruct(P);
        Pb.Tx.NTheta = num_beams(b);
        Pb.Tx.SinTheta = linspace(-Pb.Tx.SinThMax,Pb.Tx.SinThMax,Pb.Tx.NTheta);
        Pb.Tx.Theta = asin(Pb.Tx.SinTheta);

        b_BF = m_BF{b};
        b_DA = data_DA{b};
        if m == 1
            max_peak_DAS = max(db(abs(b_BF(:))));
        end

        thetaRange = Pb.Tx.Theta;
        if strcmp(methods_set(m),'IAA-MBMB-Upsampled')
            thetaRange = linspace(Pb.Tx.Theta(1), Pb.Tx.Theta(end), 500);
        end
        
        if enable_plots
            warning('off')
            [scanConvertedImage, Xs, Zs] = getScanConvertedImage(b_BF, ...
                thetaRange, 1e3 * b_DA.Radius, 2024, 2024, 'spline');
            warning('on')
            img = db(abs(scanConvertedImage));
            img = img - max_peak_DAS;

%             figure(1);
            subplot(1,2,1)
            imagesc(Xs, Zs, img)
            xlabel('azimuth [mm]');
%             imagesc(rad2deg(thetaRange),1e3 * b_DA.Radius, db(abs(b_BF)));
%             xlabel('angle [deg]');

%             xlim([-15 15])
            ylim([34 52])
%             caxis([-130  -70]);
            caxis([-50  0]);
            colorbar
            colormap(gray)
            ylabel('range [mm]');
            title(['Method :', methods_set{m}, ', No Beams: ', ...
                int2str(num_beams(b))])
        end

        thetaRange = rad2deg(thetaRange);
        % beampatterns - Assumes 2 points
        bf_img = db(abs(b_BF))  - max_peak_DAS;
        separation = P.Tx.FocRad + 5  * 1e-3; % -> 45mm
        z_sep = find(b_DA.Radius >= separation, 1);
        p1_bp = max(bf_img(1:z_sep, :), [], 1);
        p2_bp = max(bf_img(z_sep+1:end, :), [], 1);

        warning('off')
        [p1db, p1deg] = findpeaks(p1_bp, thetaRange, ...
            'SortStr','descend');
        p1_3dbline = ones(size(thetaRange)) .* (p1db(1) - 3);
        [p1_x, p1_y, ~, ~] = intersections(thetaRange, p1_bp, ...
            thetaRange, p1_3dbline, 1);
        if length(p1_x) < 2
            p1_width = NaN;
        else
            p1_width = p1_x(2) - p1_x(1);
        end

        [p2db, p2deg] = findpeaks(p2_bp, thetaRange, ...
            'SortStr','descend');
        p2_3dbline = ones(size(thetaRange)) * (p2db(1) - 3);
        [p2_x, p2_y, ~, ~] = intersections(thetaRange, p2_bp, ...
            thetaRange, p2_3dbline, 1);
        if length(p2_x) < 2
            p2_width = NaN;
        else
            p2_width = p2_x(2) - p2_x(1);
        end
        warning('on')
        points_3dBwidth{m,b} = [p1_width, p2_width];
        
        if enable_plots
    %             figure(2); 
            subplot(1,2,2)
            p = plot(thetaRange, p1_bp, 'b', ...
                thetaRange, p2_bp, 'r', 'LineWidth', 2);
            linestyle_list = {'-.','--','-',':'};
    %             markers_list = {'+','x','diamond','o'};
            for pidx=1:length(p)
    %                 p(pidx).Marker = markers_list{pidx};
                p(pidx).LineStyle = linestyle_list{pidx};
            end
            p1_title = strcat('P1-', ...
                int2str(data_phantoms.positions(1,3)*1000), 'mm range');
            p2_title = strcat('P2-', ...
                int2str(data_phantoms.positions(2,3)*1000), 'mm range');

            line('XData', [p1_x(1) p1_x(2)], 'YData', [p1_y(1) p1_y(1)], ...
                'LineWidth', 2, 'LineStyle', '-', 'Color', [0,0,0]+0.4)
            p1_l = strcat('P1-3dB width= ', num2str(p1_width, 3), 'dB');
            line('XData', [p2_x(1) p2_x(2)], 'YData', [p2_y(1) p2_y(1)], ...
                'LineWidth', 2, 'LineStyle', '-', 'Color', [0,0,0]+0.2)
            p2_l = strcat('P1-3dB width= ', num2str(p2_width, 3), 'dB');

            y_min = min([p1_y(1), p2_y(1)]) - 3;
            x_min = min([p1_x(1), p2_x(1)]) - 1;
            x_max = max([p1_x(2), p2_x(2)]) + 1;
%             xlim([-5 15])
            xlim([x_min x_max])
            ylim([y_min 0])
            xlabel('angle [deg]');
            ylabel('gain [dB]');
            legend({p1_title, p2_title, p1_l, p2_l}, 'Location', 'best');
            pause
        end
    end
end