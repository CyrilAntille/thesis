%% 2.1: Motion between frames
addpath '/uio/hume/student-u04/cyrila/Documents/MATLAB/MasterThesis/Field_II_ver_3_24' -end
% addpath 'C:\Users\Cyril\Documents\MATLAB\Field_II_ver_3_24' -end

save_all_data = true; % Can take multiple GB of memory
enable_plots = false;
grayscale_plots = true;

bfm = [0,1,3,4]; % 0=DAS, 1=MV, 2=MV-Multibeam, 3=IAA, 4=IAA_500pts
methods_set = {'DAS','MV','IAA','IAA 500pts'}; % Must correspond to bfm

disable_multiprocess = false; % Must disable automatic creation of 
% parallel pool for parfor (in parallel preferences).

%% 0. Load data
% If raw/DA/BF data exists, the script won't recreate the data.
% -> Can just load data from file if available

% load 2_1_speckle.mat % Loads raw speckle data
% bf_data_files = {'2_1_61_71.mat', '2_1_81_91.mat'};
% [ shift_type, num_beams, data_phantoms, data_DA, data_BF ] ...
%     = mergeData(bf_data_files, false); % Loads and merges DA(S) (and BF) data
% shift_type.shift = 0.5;
% shift_type.num_shifts = 2;
% Note: For now shift_type properties values are not saved (MATLAB bug)

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
if exist('speckle_raw', 'var') ~= 1
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
    shift_type = ShiftType.RadialVar;
    shift_type.num_shifts = 4;
    shift_type.shift = 1/2;
    num_beams = 81:10:91;
    
    pts_theta = [0]; % Add a theta (in degrees) for each point
    pts_range = [P.Tx.FocRad]; % Add a range (in m) for each point
    scat_pts = zeros([length(pts_theta) 3]);
    for pidx = 1:length(pts_theta)
        pt = [sind(pts_theta(pidx)) 0 cosd(pts_theta(pidx))] ...
            * pts_range(pidx);
        scat_pts(pidx, :) = pt;
    end
    original_phantom = PointPhantom(scat_pts, 1+db2mag(20));
    % -> pts 20dB over speckle
    
    data_phantoms = cell([1, length(num_beams)]);
    data_DA = cell([1, length(num_beams)]);
    for b=1:length(num_beams)
        Pb = copyStruct(P);
        Pb.Tx.NTheta = num_beams(b);
        Pb.Tx.SinTheta = linspace(-Pb.Tx.SinThMax, Pb.Tx.SinThMax, Pb.Tx.NTheta);
        Pb.Tx.Theta = asin(Pb.Tx.SinTheta);
        fprintf('\nNTheta: %d.\n', Pb.Tx.NTheta)
        nstart = tic;
        
        b_phantoms = cell([1, shift_type.num_shifts]);
        b_DA = cell([1, shift_type.num_shifts]);
        shifts = shift_type.getShifts(Pb);
        parfor s=1:shift_type.num_shifts
            s_phantom = shift_type.shiftPositions(original_phantom, shifts(s));
            b_phantoms{s} = s_phantom;
            s_raw = CalcRespAll(Pb, s_phantom);
            if add_speckle
                img = speckle_raw.image;
                % Assumes scatterers in speckle -> s_raw.image smaller than img.
                img(1:size(s_raw.image, 1),:) = img(1:size(s_raw.image, 1),:) + s_raw.image;
                s_raw.image = img;
            end
            b_DA{s} = BeamformAll(Pb, s_raw);
        end
        data_DA{b} = b_DA;
        data_phantoms{b} = b_phantoms;
        nend = toc(nstart);
        fprintf('Time: %d minutes and %f seconds\b', floor(nend/60), rem(nend,60))
    end
    if save_all_data
        output_file = ['2_1_' num2str(num_beams(1)) '_' ...
            num2str(num_beams(end)) '.mat'];
        fprintf('\nNSaving DA(S) data to file.')
        save(output_file, 'shift_type', 'num_beams', 'data_phantoms', ...
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
        for b=1:length(num_beams)
            Pb = copyStruct(P);
            Pb.Tx.NTheta = num_beams(b);
            Pb.Tx.SinTheta = linspace(-Pb.Tx.SinThMax,Pb.Tx.SinThMax,Pb.Tx.NTheta);
            Pb.Tx.Theta = asin(Pb.Tx.SinTheta);
            
            b_DA = data_DA{b};
            b_BF = cell([1, shift_type.num_shifts]);
            parfor s=1:shift_type.num_shifts
                b_BF{s} = ComputeBF(b_DA{s}.image, Pb, bf_method);
            end
            m_BF{b} = b_BF;
        end
        data_BF{m} = m_BF;
        mend = toc(mstart);
        fprintf('%s: %d minutes and %f seconds\n', methods_set{m}, ...
            floor(mend/60), rem(mend,60))
    end
    if save_all_data
        output_file = ['2_1_' num2str(num_beams(1)) '_' ...
            num2str(num_beams(end)) '.mat'];
        fprintf('\nNSaving beamformed data to file.')
        save(output_file, 'data_BF', '-append')
    end
end
fprintf('\n============================================================\n')

%% Estimate scalloping loss
fprintf('Scalloping loss estimation and plotting\b')
beam_max = zeros([length(methods_set) length(num_beams)]);
beam_mean = zeros([length(methods_set) length(num_beams)]);

theta_max = asin(P.Tx.SinThMax);
for m=1:length(methods_set)
    m_BF = data_BF{m};
    fprintf('----------------------\n')
    fprintf('Beamforming method: %s.\n', methods_set{m})
    for b=1:length(num_beams)
        b_BF = m_BF{b};
        b_DA = data_DA{b};
        num_pts = size(data_phantoms{1}{1}.positions, 1);
        pts_ampl = ones(shift_type.num_shifts, num_pts) * Inf;
        for s=1:shift_type.num_shifts
            s_BF = b_BF{s};
            s_DA = b_DA{s};
            speckle_pts = data_phantoms{b}{s};
            [Theta,Phi,R] = cart2sph(speckle_pts.positions(:,3),...
                speckle_pts.positions(:,1),speckle_pts.positions(:,2));
            
            radiusRange = s_DA.Radius;
            thetaRange = linspace(-theta_max, theta_max, num_beams(b));
            if strcmp(methods_set(m),'IAA 500pts')
                thetaRange = linspace(-theta_max, theta_max, 500);
            end
            s_pts = pts_ampl(s,:);
            for p=1:length(Theta)
                pt = Theta(p);
                pr = R(p);
                s_pts(p) = Inf;
                if pt < thetaRange(1) || pt > thetaRange(end) || ...
                        pr < radiusRange(1) || pr > radiusRange(end)
                    continue;
                end
                [~, tidx] = min(abs(pt-thetaRange));
                [~, ridx] = min(abs(pr-radiusRange));
                s_pts(p) = db(abs(s_BF(ridx,tidx)));
            end
            pts_ampl(s,:) = s_pts;
        end
        max_loss = 0; mean_loss = 0; num_pts = 0;
        for p=1:size(pts_ampl, 2)
            if max(pts_ampl(:,p)) == Inf
                continue % -> point not in beamformed area
            end
            mloss = max(pts_ampl(:,p)) - min(pts_ampl(:,p));
            max_loss = max([max_loss, mloss]);
            mean_loss = mean_loss + mloss;
            num_pts = num_pts + 1;
        end
        beam_max(m,b) = max_loss;
        beam_mean(m,b) = mean_loss/num_pts;
    end
end

%% Plots
if enable_plots
    % Max loss
    figure;
    p = plot(num_beams, beam_max', 'LineWidth', 2);
    linestyle_list = {':','-','--','-.'};
    markers_list = {'+','x','diamond','o'};
    grayscale_alpha = linspace(0,0.4,length(p));
    for pidx=1:length(p)
        p(pidx).Marker = markers_list{pidx};
        p(pidx).LineStyle = linestyle_list{pidx};
        if grayscale_plots
            p(pidx).Color = [0,0,0]+grayscale_alpha(pidx);
        end
    end

    line_clr = 'green';
    if grayscale_plots
        line_clr = [0,0,0]+0.6;
    end
    line('XData', [num_beams(1) num_beams(end)], 'YData', [1 1], ...
        'LineWidth', 2, 'LineStyle', '-.', 'Color', line_clr);
    legend([methods_set '1dB threshold'], 'Location', 'best');
    ylabel('Max scalloping loss [dB]');
    xlabel('Number of transmitted beams');
    if length(num_beams) > 1
        xlim([num_beams(1) num_beams(end)]) 
    end
    % Mean loss
    figure;
    p = plot(num_beams, beam_mean', 'LineWidth', 2);
    linestyle_list = {':','-','--','-.'};
    markers_list = {'+','x','diamond','o'};
    grayscale_alpha = linspace(0,0.4,length(p));
    for pidx=1:length(p)
        p(pidx).Marker = markers_list{pidx};
        p(pidx).LineStyle = linestyle_list{pidx};
        if grayscale_plots
            p(pidx).Color = [0,0,0]+grayscale_alpha(pidx);
        end
    end

    line_clr = 'green';
    if grayscale_plots
        line_clr = [0,0,0]+0.6;
    end
    l = line('XData', [num_beams(1) num_beams(end)], 'YData', [1 1], ...
        'LineWidth', 2, 'LineStyle', '-.', 'Color', line_clr);
    legend([methods_set '1dB threshold'], 'Location', 'best');
    ylabel('Mean scalloping loss [dB]');
    xlabel('Number of transmitted beams');
    if length(num_beams) > 1
        xlim([num_beams(1) num_beams(end)]) 
    end

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
            for s=1:shift_type.num_shifts
                s_BF = b_BF{s};
                s_DA = b_DA{s};

                thetaRange = Pb.Tx.Theta;
                if strcmp(methods_set(m),'IAA 500pts')
                    thetaRange = linspace(Pb.Tx.Theta(1), Pb.Tx.Theta(end), 500);
                end
                warning('off')
                [scanConvertedImage, Xs, Zs] = getScanConvertedImage(s_BF, ...
                    thetaRange, 1e3 * s_DA.Radius, 2024, 2024, 'spline');
                warning('on')
                img = scanConvertedImage./max(scanConvertedImage(:));
                img = db(abs(img));
                imagesc(Xs, Zs, img)
                caxis([-60 0]) %TODO: See if remove normalization of image
                xlabel('azimuth [mm]');

%                 imagesc(rad2deg(thetaRange),1e3 * s_DA.Radius, db(abs(s_BF)));
%                 caxis([-130  -70]);
%                 xlim([-15 15])
%                 ylim([35 45])
%                 xlabel('angle [deg]');

                colorbar
                if grayscale_plots
                    colormap(gray)
                else 
                    colormap(jet)
                end 
                ylabel('range [mm]');
                title(['Method :', methods_set{m}, ', No Beams: ', ...
                    int2str(num_beams(b)), ', Shift No: ', int2str(s-1)])
                pause
            end
        end
    end
end
