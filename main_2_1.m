%% 2.1: Motion between frames
addpath '/uio/hume/student-u04/cyrila/Documents/MATLAB/MasterThesis/Field_II_ver_3_24' -end

save_all_data = true; %Can take multiple GB of memory
grayscale_plots = true;

bfm = [0,1,3,4]; % 0=DAS, 1=MV, 2=MV-Multibeam, 3=IAA, 4=IAA_500pts
methods_set = {'DAS','MV','IAA','IAA 500pts'}; %Must correspond to bfm

disable_multiprocess = false; % Must disable automatic creation of 
% parallel pool for parfor (in parallel preferences).
if disable_multiprocess
    pool = struct; pool.NumWorkers = 1;
else
    if isempty(gcp('nocreate'))
        parpool;
    end
    pool = gcp;
end

%% 1. Create raw data
fprintf('Initializing Field II for all workers.\n')
field_init(0);
parfor w=1:pool.NumWorkers
    field_init(0);
end
fprintf('\n============================================================\n')

create_speckle = true;
if create_speckle && exist('speckle_raw', 'var') ~= 1
    save_speckle = true;
    fprintf('Creating raw speckle image. This can take hours if many points.\n')
    P = Parameters();
    P.Seed = 42;
    P.NumPoints = 10^6;
    speckle_phantom = PlanewaveVesselPhantom(P,0,P.NumPoints,P.Seed);
    speckle_raw = CalcRespAll(P, speckle_phantom);
    
    if save_speckle
        fprintf('\nNSaving speckle raw data to file.')
        save -v7.3 2_1_speckle.mat P speckle_raw speckle_phantom;
    end
end
fprintf('\n============================================================\n')

add_speckle = true;
if exist('speckle_raw', 'var') ~= 1
    P = Parameters();
    add_speckle = false;
    speckle_raw = [];
end
if exist('data_DA', 'var') ~= 1
    fprintf('Creating raw and DA(S) images. This can take hours if many beam setups (NThetas).\n')
    shift_type = ShiftType.RadialVar;
    shift_type.num_shifts = 4;
    shift_type.shift = 1/2;
    num_beams = 61:10:91;
    
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
%     data_raw = cell([1, length(num_beams)]);
    data_DA = cell([1, length(num_beams)]);
    for b=1:length(num_beams)
        Pb = copyStruct(P);
        Pb.Tx.NTheta = num_beams(b);
        Pb.Tx.SinTheta = linspace(-Pb.Tx.SinThMax, Pb.Tx.SinThMax, Pb.Tx.NTheta);
        Pb.Tx.Theta = asin(Pb.Tx.SinTheta);
        fprintf('\nNTheta: %d.\n', Pb.Tx.NTheta)
        nstart = tic;
        
        b_phantoms = cell([1, shift_type.num_shifts]);
%         b_raw = cell([1, shift_type.num_shifts]);
        b_DA = cell([1, shift_type.num_shifts]);
        shifts = shift_type.getShifts(Pb);
        parfor s=1:shift_type.num_shifts
            s_phantom = shift_type.shiftPositions(original_phantom, shifts(s));
            b_phantoms{s} = s_phantom;
            s_raw = CalcRespAll(Pb, s_phantom);
            if add_speckle
                if ~isequal(size(speckle_raw.image), size(s_raw.image))
                    error('Speckle image size differs from points image size! Exiting script.')
                end
                s_raw.image = s_raw.image + speckle_raw.image;
            end
%             b_raw{s} = s_raw;
            b_DA{s} = BeamformAll(Pb, s_raw);
        end
        data_DA{b} = b_DA;
        data_phantoms{b} = b_phantoms;
        nend = toc(nstart);
        fprintf('Time: %d minutes and %f seconds\b', floor(nend/60), rem(nend,60))
    end
    if save_all_data
        fprintf('\nNSaving DA(S) data to file.')
        save -v7.3 2_1_all_data.mat shift_type num_beams data_phantoms data_DA;
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
        fprintf('\nNSaving beamformed data to file.')
        save -v7.3 -append 2_1_all_data.mat data_BF;
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
        num_pts = size(data_phantoms{1}.positions, 1);
        pts_ampl = ones(shift_type.num_shifts, num_pts) * Inf;
        parfor s=1:shift_type.num_shifts
            s_BF = b_BF{s};
            s_DA = b_DA{s};
            speckle_pts = data_phantoms{s};
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
%                 if pt < thetaRange(1) || pt > thetaRange(end) || ...
                if pt < deg2rad(-14) || pt > deg2rad(14) || ...
                        pr < 36/1000 || pr > 74/1000
%                         pr < radiusRange(1) || pr > radiusRange(end)
                    continue;
                end
                [~, tidx] = min(abs(pt-thetaRange));
                [~, ridx] = min(abs(pr-radiusRange));
                s_pts(p) = s_BF(ridx,tidx);
%                 if s_pts(p) < -100
%                     s_BF(ridx-2:ridx+2,tidx-2:tidx+2) = - Inf;
%                 end
            end
            pts_ampl(s,:) = s_pts;
%             n_BF{s} = s_BF;
        end
%         data_BF{m}{b} = n_BF;
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

%% Max loss
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
xlim([num_beams(1) num_beams(end)]) 

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
xlim([num_beams(1) num_beams(end)])

%% Plot images
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
            
            [scanConvertedImage, Xs, Zs] = getScanConvertedImage(s_BF.', ...
                Pb.Tx.Theta, 1e3 * s_DA.Radius, 2024, 2024, 'spline');
%             imagesc(rad2deg(Pb.Tx.Theta),1e3 * s_DA.Radius, s_BF);
            img = scanConvertedImage./max(scanConvertedImage(:));
            img = db(abs(img));
            imagesc(Xs, Zs, img)
            
%             caxis([-130  -70]);
            caxis([-60 0])
%             xlim([-15 15])
%             ylim([35 75])
            
            colorbar
            if grayscale_plots
                colormap(gray)
%                 colormap(flipud(colormap)) 
            else 
                colormap(jet)
            end 
            ylabel('range [mm]');
            xlabel('angle [deg]');
            title(['Method :', methods_set{m}, ', No Beams: ', ...
                int2str(num_beams(b)), ', Shift No: ', int2str(s-1)])
            pause
        end
    end
end
