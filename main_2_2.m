%% 2.2: Motion between beams
addpath '/uio/hume/student-u04/cyrila/Documents/MATLAB/MasterThesis/Field_II_ver_3_24' -end
% addpath 'C:\Users\Cyril\Documents\MATLAB\Field_II_ver_3_24' -end

save_all_data = false; % Can take multiple GB of memory
enable_plots = false;
grayscale_plots = true;

bfm = [0,1,3,4]; % 0=DAS, 1=MV, 2=MV-Multibeam, 3=IAA, 4=IAA_500pts
methods_set = {'DAS','MV','IAA','IAA 500pts'}; %Must correspond to bfm

disable_multiprocess = false; % Must disable automatic creation of 
% parallel pool for parfor (in parallel preferences).

%% 0. Load data
% If raw/DA/BF data exists, the script won't recreate the data.
% -> Can just load data from file if available

% load 2_1_speckle.mat % Loads raw speckle data
% bf_data_files = {'2_1_61_91.mat', '2_1_101_131.mat'};
% [ shift_type, num_beams, data_phantoms, data_DA, data_BF ] ...
%     = mergeData(bf_data_files, false); % Loads and merges DA(S) (and BF) data

%% 1. Create Speckle raw data
field_init(0);

create_speckle = true; 
if create_speckle && exist('speckle_raw', 'var') ~= 1
    save_speckle = true;
    fprintf('Creating raw speckle image. This can take hours if many points.\n')
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
    num_beams = 61:10:71;
    
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
    
            s_phantom = shift_type.shiftPositions(original_phantom, shifts(s));
            b_phantoms{s} = s_phantom;
            s_raw = CalcRespAll(Pb, s_phantom);
            if add_speckle
                img = speckle_raw.image;
                % Assumes scatterers in speckle -> s_raw.image smaller than img.
                img(1:size(s_raw.image, 1),:) = img(1:size(s_raw.image, 1),:) + s_raw.image;
                s_raw.image = img;
            end
    parfor b=1:length(num_beams)
        b_DA = cell([1, num_beams(b)]);
        
        Pb = copyStruct(P);
        Pb.Tx.NTheta = num_beams(b);
        Pb.Tx.SinTheta = sin_thetas(b);
        Pb.Tx.Theta = asin(Pb.Tx.SinTheta);
        shifts = shift_type.getShifts(Pb);
        sin_thetas = linspace(-P.Tx.SinThMax,P.Tx.SinThMax, num_beams(n));
        fprintf('\nNTheta: %d.\n', P.Tx.NTheta)
        nstart = tic;
        for n=1:num_beams(n)
            Pb = copyStruct(P);
            Pb.Tx.NTheta = 1;
            Pb.Tx.SinTheta = sin_thetas(b);
            Pb.Tx.Theta = asin(Pb.Tx.SinTheta);
            
            s_phantom = shift_type.shiftPositions(original_phantom, shifts(n));
            s = ceil(b / num_beams(n) * num_shifts);
            b_DA{b} = BeamformAll(Pb, data_raw{s});
        end
        n_DA = b_DA{1};
        for b=1:num_beams(n)
            radius_length = length(b_DA{b}.Radius); % This number may vary
            if length(n_DA.Radius) > radius_length
                n_DA.image = n_DA.image(:,1:radius_length,:);
                n_DA.Radius = b_DA{b}.Radius;
            elseif radius_length > length(n_DA.Radius)
                b_DA{b}.image = b_DA{b}.image(:,1:length(n_DA.Radius),:);
            end
            n_DA.image = vertcat(n_DA.image, b_DA{b}.image);
        end
        data_DA{n} = n_DA;
    end
    if save_data
        fprintf('\nNSaving DA data to file. ')
        save -v7.3 -append 2_2_all_data.mat num_beams data_DA;
    end
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
%% 2.2: Motion between beams
addpath '/uio/hume/student-u04/cyrila/Documents/MATLAB/MasterThesis/Field_II_ver_3_24' -end


save_raw_data = false;
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
    
    create_speckle = false;
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
        P.Seed = 2;
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
    if save_data
        fprintf('\nNSaving raw data to file.')
        save -v7.3 2_2_all_data.mat P img_shifts_rad img_shifts_mm data_raw data_phantoms;
    end
else
    num_shifts = length(img_shifts_rad);
    if any(img_shifts_mm)
        num_shifts = length(img_shifts_mm);
    end
end
fprintf('\n============================================================\n')

%% 2. Create DA(S) data
% Considers the scene to move of num_shifts between first and last beam.
if exist('data_DA', 'var') ~= 1
    fprintf('Creating DA(S) images. \n')
    P_shift = Parameters(180);
    num_beams = 61:10:91;
    data_DA = cell([1, length(num_beams)]);
    parfor n=1:length(num_beams)
        b_DA = cell([1, num_beams(n)]);
        sin_thetas = linspace(-P.Tx.SinThMax,P.Tx.SinThMax, num_beams(n));
        fprintf('\nNTheta: %d.\n', Pb.Tx.NTheta)
        nstart = tic;
        for b=1:num_beams(n)
            Pb = copyStruct(P);
            Pb.Tx.NTheta = 1;
            Pb.Tx.SinTheta = sin_thetas(b);
            Pb.Tx.Theta = asin(Pb.Tx.SinTheta);
            s = ceil(b / num_beams(n) * num_shifts);
            b_DA{b} = BeamformAll(Pb, data_raw{s});
        end
        n_DA = b_DA{1};
        for b=1:num_beams(n)
            radius_length = length(b_DA{b}.Radius); % This number may vary
            if length(n_DA.Radius) > radius_length
                n_DA.image = n_DA.image(:,1:radius_length,:);
                n_DA.Radius = b_DA{b}.Radius;
            elseif radius_length > length(n_DA.Radius)
                b_DA{b}.image = b_DA{b}.image(:,1:length(n_DA.Radius),:);
            end
            n_DA.image = vertcat(n_DA.image, b_DA{b}.image);
        end
        data_DA{n} = n_DA;
    end
    if save_data
        fprintf('\nNSaving DA data to file. ')
        save -v7.3 -append 2_2_all_data.mat num_beams data_DA;
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
        parfor n=1:length(num_beams)
            Pb = copyStruct(P);
            Pb.Tx.NTheta = num_beams(n);
            Pb.Tx.SinTheta = linspace(-Pb.Tx.SinThMax,Pb.Tx.SinThMax,Pb.Tx.NTheta);
            Pb.Tx.Theta = asin(Pb.Tx.SinTheta);
            n_BF = ComputeBF(data_DA{n}.image, Pb, bf_method);
        end
        data_BF{m} = m_BF;
        mend = toc(mstart);
        fprintf('%s: %d minutes and %f seconds\n', methods_set{m}, ...
            floor(mend/60), rem(mend,60))
    end
    if save_data
        fprintf('\nNSaving beamformed data to file.')
        save -v7.3 -append 2_2_all_data.mat data_BF;
    end
end
fprintf('\n============================================================\n')

%% Plot images
figure;
for m=1:length(methods_set)
    m_BF = data_BF{m};
    for n=1:length(num_beams)
        Pb = copyStruct(P);
        Pb.Tx.NTheta = num_beams(n);
        Pb.Tx.SinTheta = linspace(-Pb.Tx.SinThMax,Pb.Tx.SinThMax,Pb.Tx.NTheta);
        Pb.Tx.Theta = asin(Pb.Tx.SinTheta);

        imagesc(rad2deg(Pb.Tx.Theta),1e3 * data_DA{n}.Radius, m_BF{n});

        caxis([-130  -70]);
        xlim([-15 15])
        ylim([35 75])

        colorbar
        if grayscale_plots
            colormap(gray)
            colormap(flipud(colormap)) 
        else 
            colormap(jet)
        end 
        ylabel('range [mm]');
        xlabel('angle [deg]');
        title(['Method :', methods_set{m}, ', No Beams: ', ...
            int2str(num_beams(n))])
        pause
    end
end

%% Estimate scalloping loss - TO RETHINK
fprintf('Scalloping loss estimation and plotting\n')
beam_max = zeros([length(methods_set) length(num_beams)]);
beam_mean = zeros([length(methods_set) length(num_beams)]);

theta_max = asin(P.Tx.SinThMax);
for m=1:length(methods_set)
    m_BF = data_BF{m};
    fprintf('----------------------\n')
    fprintf('Beamforming method: %s.\n', methods_set{m})
    parfor n=1:length(num_beams)
        n_BF = m_BF{n};
        b_DA = data_DA{n};
        num_pts = length(data_phantoms{1}.positions);
        pts_ampl = ones(num_shifts, num_pts) * Inf;
        for s=1:num_shifts
            s_BF = n_BF{s};
            s_DA = b_DA{s};
            speckle_pts = data_phantoms{s};
            [Theta,Phi,R] = cart2sph(speckle_pts.positions(:,3),...
                speckle_pts.positions(:,1),speckle_pts.positions(:,2));
            
            radiusRange = s_DA.Radius;
            thetaRange = linspace(-theta_max, theta_max, num_beams(n));
            if strcmp(methods_set(m),'IAA 500pts')
                thetaRange = linspace(-theta_max, theta_max, 500);
            end
            for p=1:length(Theta)
                pt = Theta(p);
                pr = R(p);
                pts_ampl(s,p) = Inf;
%                 if pt < thetaRange(1) || pt > thetaRange(end) || ...
                if pt < deg2rad(-14) || pt > deg2rad(14) || ...
                        pr < 36/1000 || pr > 74/1000
%                         pr < radiusRange(1) || pr > radiusRange(end)
                    continue;
                end
                [~, tidx] = min(abs(pt-thetaRange));
                [~, ridx] = min(abs(pr-radiusRange));
                pts_ampl(s,p) = s_BF(ridx,tidx);
            end
        end
        max_loss = 0; mean_loss = 0; num_pts = 0;
        for p=1:length(pts_ampl)
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
legend(methods_set, '1dB threshold');
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
legend(methods_set, '1dB threshold');
ylabel('Mean scalloping loss [dB]');
xlabel('Number of transmitted beams');
xlim([num_beams(1) num_beams(end)])

