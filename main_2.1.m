%% 2.1: Motion between frames
addpath '/uio/hume/student-u04/cyrila/Documents/MATLAB/MasterThesis/Field_II_ver_3_24' -end

save_raw_data = true;
save_all_data = save_raw_data && false; %Can take multiple GB of memory
grayscale_plots = false;

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
parfor w=1:pool.NumWorkers
    field_init(0);
end

fprintf('\n============================================================\n')
if exist('data_raw', 'var') ~= 1
    fprintf('Creating raw images. This can take many hours if many points.\n')
    P = Parameters();
    
    create_speckle = true;
    curved_shift = true; % if true, uses img_shifts_rad, else uses img_shifts_mm
    num_shifts = 4;
    
    img_shifts_rad = zeros([1 num_shifts]);
    img_shifts_mm = zeros([1 num_shifts]);
    if curved_shift
        P_shift = Parameters(181);
        shft = (P_shift.Tx.Theta(2) - P_shift.Tx.Theta(1)) / 2; 
        img_shifts_rad = (0:num_shifts-1) * shft;
    else
        % Note: sin(1 degree) * 40mm = about 0.7mm
        % Theta range [-17.5, 17.5] degrees -> for 81 beams: beams distance
        % is 35 / 81 = 0.432 degrees
        shft = 3 / 1e4;
        img_shifts_mm = (0:num_shifts-1) * shft;
        [~,img_shifts_mm,~] = sph2cart(img_shifts_rad, 0, P.Tx.FocRad);
    end
    
    data_raw = cell([1, num_shifts]);
    data_phantoms = cell([1, num_shifts]);
    if create_speckle
        P.Seed = 42;
        P.NumPoints = 10^6;
        original_phantom = PlanewaveVesselPhantom(P,0,P.NumPoints,P.Seed);
    else % Scatterer points
        pts_theta = [0]; % Add a theta (in degrees) for each point
        pts_range = [P.Tx.FocRad]; % Add a range (in m) for each point
        scat_pts = zeros([length(pts_theta) 3]);
        for pidx = 1:length(pts_theta)
            pt = [sind(pts_theta(pidx)) 0 cosd(pts_theta(pidx))] ...
                * pts_range(pidx);
            scat_pts(pidx, :) = pt;
        end
        original_phantom = PointPhantom(scat_pts);
    end
    parfor s=1:num_shifts
        s_phantom = copyStruct(original_phantom);
        if curved_shift
            [Theta,Phi,R] = cart2sph(s_phantom.positions(:,3),...
                s_phantom.positions(:,1),s_phantom.positions(:,2));
            Theta = Theta + img_shifts_rad(s);
            [z,x,y] = sph2cart(Theta,Phi,R);
            s_phantom.positions = [x,y,z];
        else %linear shift (along x)
            s_phantom.positions = [s_phantom.positions(:,1) + img_shifts_mm(s),...
                s_phantom.positions(:,2), s_phantom.positions(:,3)];
        end
        data_raw{s} = CalcRespAll(P,s_phantom);
        data_phantoms{s} = s_phantom;
    end
    if save_raw_data
        fprintf('\nNSaving raw data to file.')
        save -v7.3 2_1_all_data.mat P img_shifts_rad img_shifts_mm data_raw data_phantoms;
    end
else
    num_shifts = length(img_shifts_rad);
    if any(img_shifts_mm)
        num_shifts = length(img_shifts_mm);
    end
end
fprintf('\n============================================================\n')

%% 2. Create DA(S) data
if exist('data_DA', 'var') ~= 1
    fprintf('Creating DA(S) images. This can take hours if many beam setups (NThetas).\n')
    num_beams = 61:10:91;
    data_DA = cell([1, length(num_beams)]);
    for n=1:length(num_beams)
        n_DA = cell([1, num_shifts]);
        Pn = copyStruct(P);
        Pn.Tx.NTheta = num_beams(n);
        Pn.Tx.SinTheta = linspace(-Pn.Tx.SinThMax,Pn.Tx.SinThMax,Pn.Tx.NTheta);
        Pn.Tx.Theta = asin(Pn.Tx.SinTheta);
        fprintf('\nNTheta: %d.\n', Pn.Tx.NTheta)
        nstart = tic;
        parfor s=1:num_shifts
            n_DA{s} = BeamformAll(Pn,data_raw{s});
        end
        data_DA{n} = n_DA;
        nend = toc(nstart);
        fprintf('Time: %d minutes and %f seconds\n', floor(nend/60), rem(nend,60))
    end
    if save_all_data
        fprintf('\nNSaving DA data to file. ')
        save -v7.3 -append 2_1_all_data.mat num_beams data_DA;
    end
end
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
        for n=1:length(num_beams)
            Pn = copyStruct(P);
            Pn.Tx.NTheta = num_beams(n);
            Pn.Tx.SinTheta = linspace(-Pn.Tx.SinThMax,Pn.Tx.SinThMax,Pn.Tx.NTheta);
            Pn.Tx.Theta = asin(Pn.Tx.SinTheta);
            
            n_DA = data_DA{n};
            n_BF = cell([1, num_shifts]);
            parfor s=1:num_shifts
                n_BF{s} = ComputeBF(n_DA{s}.image, Pn, bf_method);
            end
            m_BF{n} = n_BF;
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
fprintf('Scalloping loss estimation and plotting\n')
beam_max = zeros([length(methods_set) length(num_beams)]);
beam_mean = zeros([length(methods_set) length(num_beams)]);

theta_max = asin(P.Tx.SinThMax);
for m=1:length(methods_set)
    m_BF = data_BF{m};
    fprintf('----------------------\n')
    fprintf('Beamforming method: %s.\n', methods_set{m})
    for n=1:length(num_beams)
        n_BF = m_BF{n};
        n_DA = data_DA{n};
        num_pts = size(data_phantoms{1}.positions, 1);
        pts_ampl = ones(num_shifts, num_pts) * Inf;
        parfor s=1:num_shifts
            s_BF = n_BF{s};
            s_DA = n_DA{s};
            speckle_pts = data_phantoms{s};
            [Theta,Phi,R] = cart2sph(speckle_pts.positions(:,3),...
                speckle_pts.positions(:,1),speckle_pts.positions(:,2));
            
            radiusRange = s_DA.Radius;
            thetaRange = linspace(-theta_max, theta_max, num_beams(n));
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
%         data_BF{m}{n} = n_BF;
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
        beam_max(m,n) = max_loss;
        beam_mean(m,n) = mean_loss/num_pts;
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
    for n=1:length(num_beams)
        Pn = copyStruct(P);
        Pn.Tx.NTheta = num_beams(n);
        Pn.Tx.SinTheta = linspace(-Pn.Tx.SinThMax,Pn.Tx.SinThMax,Pn.Tx.NTheta);
        Pn.Tx.Theta = asin(Pn.Tx.SinTheta);

        n_BF = m_BF{n};
        n_DA = data_DA{n};
        for s=1:num_shifts
            s_BF = n_BF{s};
            s_DA = n_DA{s};
            
            
            
            [scanConvertedImage, Xs, Zs] = getScanConvertedImage(s_BF.', ...
                Pn.Tx.Theta, 1e3 * s_DA.Radius, 2024, 2024);
%             imagesc(rad2deg(Pn.Tx.Theta),1e3 * s_DA.Radius, s_BF);
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
                int2str(num_beams(n)), ', Shift No: ', int2str(s-1)])
            pause
        end
    end
end
