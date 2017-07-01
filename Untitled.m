nsds = [0, 10, 96];
mainP = MainParameters();
mainP.speckle_load = false;
mainP.normalize_bfim = false;
mainP.save_plots = true;

for n=1:length(nsds)
    mainP.nsd = nsds(n);
    mainP = mainP.createOutputDir();
    main_init
    plotBFImages(mainP, data_DA, data_BF)
    clearvars -except mainP nsds
end

%%

for m=1:length(data_BF)
    m_BF = data_BF{m};
    for s=1:length(m_BF)
        s_BF = m_BF{s};
        s_DA = data_DA{s};
        max_sBF = max(s_BF(:));
        sidx = find(s_BF == max_sBF, 1);
        srad = s_DA.Radius(mod(sidx, size(s_BF, 1)));
        stheta = mainP.P.Tx.Theta(floor(sidx / size(s_BF, 1)) + 1);
        fprintf('\n%s: shift %d, rad: %0.4f, theta: %0.3f, val: %0.3f\n', ...
            mainP.methods_set{m}, s, srad, rad2deg(stheta), max_sBF)
    end
end

%% Used for 2_1_3 plots to compare speckle and non-speckle scenarios
source_folder = '../images/main_2_1_april/2_1_3_final/';
scenario_folder = {'nospeckle/', 'speckle_2/', 'speckle_42/'};
scenario_name = {'nospeckle', 'seed-2', 'seed-42'};
result_file = 'loss_beams_p40_results_2_1_3.mat';

% 1. Load and modify mainP
load(strcat(source_folder, scenario_folder{1}, 'mainP.mat'));
mainP = obj;
methods_set = mainP.methods_set;
methods_mainP = cell([1, length(methods_set)]);
for m=1:length(methods_set)
    mainP.save_plots = true;
    mainP.save_folder = strcat('../data/test/', methods_set{m}, '/');
    mainP.methods_set = cell([1 length(scenario_name)]);
    for s=1:length(scenario_name)
        mainP.methods_set{s} = strcat(methods_set{m}, '-', scenario_name{s});
    end
    mainP = mainP.createOutputDir();
    methods_mainP{m} = mainP;
end

% 2. Load and modify pts_gain
ext_folder = {'nospeckle_ext_final/', 'speckle_2_ext_fake/', 'speckle_42_ext/'};
ext = 8;

load(strcat(source_folder, scenario_folder{1}, result_file));
[p,~,b,s] = size(pts_gain);
gain_size = [p,length(scenario_folder),b+ext,s];
methods_gain = cell([1, length(methods_set)]);
for m=1:length(methods_set)
    methods_gain{m} = zeros(gain_size);
end
for f=1:length(scenario_folder)
    load(strcat(source_folder, scenario_folder{f}, result_file));
    for m=1:length(methods_set)
%         methods_gain{m}(:,f,:,:) = pts_gain(:,m,:,:);
        methods_gain{m}(:,f,1:b,:) = pts_gain(:,m,:,:);
    end
    load(strcat(source_folder, ext_folder{f}, result_file));
    for m=1:length(methods_set)
        methods_gain{m}(:,f,b+1:b+ext,:) = pts_gain(:,m,2:end,:);
    end
end
num_beams = cat(2, 11:10:151, 251:100:951);

%% main_2_1_3 plots
for mp=1:length(methods_set)
    mainP = methods_mainP{mp};
    pts_gain = methods_gain{mp};
    
    % Plots
    linestyle_list = {'-.','--','-',':'};
    markers_list = {'+','x','diamond','o'};
    colors_list = {'b','r','g','k','m','c'};
    if mainP.save_plots
        figure('units','normalized','position',[.2 .3 .5 .3],'Visible','off')
    else
        figure;
    end

    for p=1:length(mainP.pts_range)
        clf();
        scallop_loss = zeros(length(mainP.methods_set), length(num_beams));
        for m=1:length(mainP.methods_set)
            for b=1:length(num_beams)
                scallop_loss(m,b) = max(pts_gain(p, m, b, :)) - ...
                    min(pts_gain(p, m, b, :));
            end
        end
        pl = plot(num_beams, scallop_loss', 'LineWidth', 2);
        for pidx=1:length(pl)
            pl(pidx).Marker = markers_list{pidx};
            pl(pidx).LineStyle = linestyle_list{pidx};
            pl(pidx).Color = colors_list{pidx};
        end
        hold on;
        clr_idx = length(pl) + 1;
        line('XData', [num_beams(1) num_beams(end)], 'YData', [1 1], ...
            'LineWidth', 2, 'LineStyle', linestyle_list{mod(clr_idx,length(linestyle_list))+1}, ...
            'Color', colors_list{mod(clr_idx,length(colors_list))+1});
        hold off;
        legend([mainP.methods_set '1dB threshold'], 'Location', 'best');
        ylabel('Max scalloping loss [dB]');
        xlabel('Number of transmitted beams');
        if length(num_beams) > 1
            xlim([num_beams(1) num_beams(end)]) 
        end
    %     ylim([0 4])
        if mainP.save_plots
            mainP.files_prefix = strcat('loss_beams_p', ...
                int2str(mainP.pts_range(p)), '_');
            saveas(gcf, mainP.outputFileName('png'), 'png')
            saveas(gcf, mainP.outputFileName('fig'), 'fig')
            save(strcat(mainP.save_folder, mainP.files_prefix, ...
                'results_2_1_3.mat'), 'num_beams', 'pts_gain', '-v7.3')
        else
            pause
        end
    end
    mainP.files_prefix = '';
    close

    % Plots fit
    linestyle_list = {'-.','--','-',':'};
    markers_list = {'+','x','diamond','o'};
    colors_list = {'b','r','g','k','m','c'};
    if mainP.save_plots
        figure('units','normalized','position',[.2 .3 .5 .3],'Visible','off')
    else
        figure;
    end

    for p=1:length(mainP.pts_range)
        clf();
        hold on;
        scallop_loss = zeros(length(mainP.methods_set), length(num_beams));
        loss_fit = cell([1 length(mainP.methods_set)]);
        for m=1:length(mainP.methods_set)
            for b=1:length(num_beams)
                scallop_loss(m,b) = max(pts_gain(p, m, b, :)) - ...
                    min(pts_gain(p, m, b, :));
            end
            loss_fit{m} = fit(num_beams', scallop_loss(m,:)', 'power2');
            scatter(num_beams, scallop_loss(m,:), colors_list{m}, ...
                'Marker',  markers_list{m}, 'LineWidth', 2);
        end

        clr_idx = length(mainP.methods_set) + 1;
        line('XData', [num_beams(1) num_beams(end)], 'YData', [1 1], ...
            'LineWidth', 2, 'LineStyle', linestyle_list{mod(clr_idx,length(linestyle_list))+1}, ...
            'Color', colors_list{mod(clr_idx,length(colors_list))+1});
        if length(num_beams) > 1
            xlim([num_beams(1) num_beams(end)]) 
        end

        for m=1:length(mainP.methods_set)
            pl = plot(loss_fit{m}, 'predobs');
            for pidx=1:length(pl)
                pl(pidx).LineStyle = linestyle_list{m};
                pl(pidx).Color = colors_list{m};
                pl(pidx).LineWidth = 2;
            end
    %         loss_legend{2*m-1} = mainP.methods_set{m};
    %         loss_legend{2*m} = strcat(mainP.methods_set{m}, '- Power fit');
    %         loss_legend{2*m} = strcat(mainP.methods_set{m}, ...
    %             '- Power fit: ', num2str(loss_fit{m}.a), '*x^', ...
    %             num2str(loss_fit{m}.b), ' + ', num2str(loss_fit{m}.c));
        end

        hold off;
        legend([mainP.methods_set '1dB threshold'], 'Location', 'best');
        ylabel('Max scalloping loss [dB]');
        xlabel('Number of transmitted beams');
        if mainP.save_plots
            mainP.files_prefix = strcat('loss_beams_fit_p', ...
                int2str(mainP.pts_range(p)), '_');
            saveas(gcf, mainP.outputFileName('png'), 'png')
            saveas(gcf, mainP.outputFileName('fig'), 'fig')
            save(strcat(mainP.save_folder, mainP.files_prefix, ...
                'results_2_1_3.mat'), 'num_beams', 'pts_gain', 'loss_fit', '-v7.3')
        else
            pause
        end
    end
    mainP.files_prefix = '';
    close
end
fprintf('Main_2_1_3 finished!')

%% Simulated probe bandwidth
P = Parameters();
NFFT = 4096;
ft = fft(P.Tx.ImpulseResp, NFFT);
db_onesided = db(abs(ft(1:NFFT/2+1)));
f = P.fs/2 * linspace(0,1,NFFT/2+1);

max_imp = max(db_onesided(:));
max_idx = find(db_onesided > max_imp - 6);

bdw = [f(max_idx(1)), f(max_idx(end))]

figure; plot(f, db_onesided)

%%
% k = 2pi/lambdac sin(theta)
% TODO: Explain optimal number of subdimensions based on image sector
% nsd = 2 * sin(theta) * N_elements
% 
% k+ = (K-1)/2 = x / nel = sin(theta)/2
% 
% theta = 17.5
% -> find x and compare to nsd/2.
% 
% plot beamspace projections (See INF5410) with corresponding image sector
% 
% k = 4*pi * 0.3

%% heavy computation
% clear all
% main_2_1_3
% clearvars -except mainP
% mainP.speckle_file = '..\data\2_1_speckle_2_10-6.mat';
% main_2_1_3
clear all
% num_beams = [65, 141, 145, 951];
num_beams = [65];
% num_beams_speckle = [75, 91, 101, 301];

mainP = MainParameters();
mainP.pts_range = [40];
mainP.pts_azimuth = [0];
mainP.NoMLA = 1;
mainP.shift_per_beam = true;
mainP.speckle_load = false;
% mainP.speckle_file = '..\data\2_1_speckle_42_10-6.mat';
% mainP.methods_set = {'DAS', 'IAA-MBMB', 'IAA-MBMB-2', 'IAA-MBMB-4'};
mainP.save_plots = true;
mainP = mainP.createOutputDir();
origin_folder = mainP.save_folder;
fprintf('\n--------------------------------------------------------------------------------------------------------------------------------------\n')

% 2.1.1: Single point linear velocity
for b=1:length(num_beams)
    mainP.num_beams = num_beams(b);
    mainP.P = mainP.copyP(mainP.num_beams);
    mainP.shift = Shift(ShiftType.LinearSpeed, 0,  mainP.num_beams, 0, 1);
    mainP.save_folder = strcat(origin_folder, '2.1.1/');
    mainP = mainP.createOutputDir();
    stats_2_2
    clearvars -except mainP num_beams origin_folder b
end
fprintf('\n-------------------------------------------- 2.1.1 finished! ----------------------------------------------------------------\n')

% 2.1.2: Single point various motion directions
directions = [-45, 45, 90];
for b=1:length(num_beams)
    for d=1:length(directions)
        mainP.num_beams = num_beams(b);
        mainP.P = mainP.copyP(mainP.num_beams);
        mainP.shift = Shift(ShiftType.LinearSpeed, 0,  mainP.num_beams, ...
            directions(d), 1);
        mainP.save_folder = strcat(origin_folder, '2.1.2/');
        mainP = mainP.createOutputDir();
        stats_2_2
        clearvars -except mainP num_beams directions origin_folder b d
    end
end
fprintf('\n-------------------------------------------- 2.1.2 finished! ----------------------------------------------------------------\n')

% 2.2.1: Two points linear velocity
% mainP.speckle_load = true;
mainP.pts_range = [40, 40];
img_azimuth_length = 2 * tan(asin(0.3)) * mainP.pts_range(1);
mainP.pts_azimuth = [0, 2.5*img_azimuth_length/65];
for b=1:length(num_beams)
    mainP.num_beams = num_beams(b);
    mainP.P = mainP.copyP(mainP.num_beams);
    mainP.shift = Shift(ShiftType.LinearSpeed, 0,  mainP.num_beams, 0, 1);
    mainP.save_folder = strcat(origin_folder, '2.2.1/');
    mainP = mainP.createOutputDir();
    stats_2_2
    clearvars -except mainP num_beams origin_folder b
end
fprintf('\n-------------------------------------------- 2.2.1 finished! ----------------------------------------------------------------\n')

% 2.2.2: Two points various motion directions
directions = [-45, 45, 90];
% directions = [45, 90];
for b=1:length(num_beams)
    for d=1:length(directions)
        mainP.num_beams = num_beams(b);
        mainP.P = mainP.copyP(mainP.num_beams);
        mainP.shift = Shift(ShiftType.LinearSpeed, 0,  mainP.num_beams, ...
            directions(d), 1);
        mainP.save_folder = strcat(origin_folder, '2.2.2/');
        mainP = mainP.createOutputDir();
        stats_2_2
        clearvars -except mainP num_beams directions origin_folder b d
    end
end
fprintf('\n-------------------------------------------- 2.2.2 finished! ----------------------------------------------------------------\n')

%% Plots merging for stats_2_2
%% Used for 2_1_3 plots to compare speckle and non-speckle scenarios
source_folder = '../output/2017-05-15_18.46.59/2_1_1/';
scenario_folder = {'2017-05-15_18.46.59/', '2017-05-15_19.34.06/', ...
    '2017-05-15_20.57.25/'};
scenario_name = {'65', '141', '145'};
result_file = 'speeds_p1_results_stats_2_2.mat';

% 1. Load and modify mainP
load(strcat(source_folder, scenario_folder{1}, 'mainP.mat'));
mainP = obj;
methods_set = mainP.methods_set;
methods_mainP = cell([1, length(methods_set)]);
for m=1:length(methods_set)
    mainP.save_plots = true;
    mainP.save_folder = strcat(source_folder, 'merged/', methods_set{m}, '/');
    mainP.methods_set = cell([1 length(scenario_name)]);
    for s=1:length(scenario_name)
        mainP.methods_set{s} = strcat(methods_set{m}, '-', scenario_name{s});
    end
    mainP = mainP.createOutputDir();
    methods_mainP{m} = mainP;
end

% 2. Load and modify pts_gain
ext_folder = {'nospeckle_ext_final/', 'speckle_2_ext_fake/', 'speckle_42_ext/'};
ext = 8;

load(strcat(source_folder, scenario_folder{1}, result_file));
[p,~,b,s] = size(pts_gain);
gain_size = [p,length(scenario_folder),b+ext,s];
methods_gain = cell([1, length(methods_set)]);
for m=1:length(methods_set)
    methods_gain{m} = zeros(gain_size);
end
for f=1:length(scenario_folder)
    load(strcat(source_folder, scenario_folder{f}, result_file));
    for m=1:length(methods_set)
%         methods_gain{m}(:,f,:,:) = pts_gain(:,m,:,:);
        methods_gain{m}(:,f,1:b,:) = pts_gain(:,m,:,:);
    end
    load(strcat(source_folder, ext_folder{f}, result_file));
    for m=1:length(methods_set)
        methods_gain{m}(:,f,b+1:b+ext,:) = pts_gain(:,m,2:end,:);
    end
end
num_beams = cat(2, 11:10:151, 251:100:951);


%% stats_2_2 num beams
mainP = MainParameters();
mainP.pts_range = [40, 40]; % Add a range (in mm) for each point
mainP.pts_azimuth = [0, 2*25.16/65];
mainP.num_beams = 65;
mainP.shift = Shift(ShiftType.LinearSpeed, 0, mainP.num_beams, 0, 1);
mainP.shift_per_beam = true;
mainP.save_plots = true;
mainP.speckle_load = false;
mainP.normalize_bfim = true;
mainP.norm_variant = 2;
mainP = mainP.createOutputDir();
orig_folder = mainP.save_folder;

num_beams = [141, 145, 951];
for b=1:length(num_beams)
    mainP.num_beams = num_beams(b);
    fprintf('moving_2_2: Running main_2_2 with %d beams.\n',num_beams(b));
%     pts_dist = 2 * 25.16 / mainP.num_beams;
%     mainP.pts_azimuth = [0, pts_dist];
    mainP.P = mainP.copyP(mainP.num_beams);
    mainP.shift = Shift(ShiftType.LinearSpeed, 0, mainP.num_beams, 0, 1);
    mainP.save_folder = orig_folder;
    mainP = mainP.createOutputDir();
    main_2_2
    clearvars -except mainP b num_beams orig_folder
end

%% stats_2_2 points dist
mainP = MainParameters();
mainP.pts_range = [40, 40]; % Add a range (in mm) for each point
mainP.num_beams = 145;
pts_dist = 2*25.16/mainP.num_beams;
mainP.pts_azimuth = [0, pts_dist];
mainP.shift = Shift(ShiftType.LinearSpeed, 0, mainP.num_beams, 0, 1);
mainP.shift_per_beam = true;
mainP.save_plots = true;
mainP.speckle_load = false;
mainP.normalize_bfim = true;
mainP.norm_variant = 2;
mainP.P = mainP.copyP(mainP.num_beams);
mainP = mainP.createOutputDir();
orig_folder = mainP.save_folder;

speeds = -0.5:1/8:0.5;
for s=1:length(speeds)
    fprintf('Stats_2_2: Running main_2_2 with speed value: %0.2f.\n', speeds(s));
    pts_shift = speeds(s) * 0.4; % 0.4 ms = 1/2 + 1 + 1/2 beam acqu. time
    diff_beams = 2;
    while abs(pts_shift) >= pts_dist / 2 % -> different 'main' beam
        diff_beams = diff_beams + sign(pts_shift);
        pts_shift = speeds(s) * diff_beams;
    end
    mainP.pts_azimuth = [0, pts_dist + pts_shift];
    fprintf('Diff beams: %d, pts dist: %0.3f mm\n', diff_beams, pts_dist + pts_shift)
    mainP.save_folder = orig_folder;
    mainP = mainP.createOutputDir();
    main_2_2
    clearvars -except mainP b speeds orig_folder pts_dist
end

%% Result - first illustration
mainP = MainParameters();
mainP.pts_range = [40, 50];
mainP.pts_azimuth = [0, 0];
mainP.num_beams = 15;
mainP.shift = Shift(ShiftType.LinearSpeed, 0, mainP.num_beams, 0, 1);
mainP.shift_per_beam = true;
mainP.speckle_load = true;
mainP.pts_gain = 0;

mainP.save_plots = false;
mainP.P = mainP.copyP(mainP.num_beams);
mainP = mainP.createOutputDir();

main_init
mainP.shift = Shift(ShiftType.RadialVar, 1/8, 17, 0, 1);
plotBFImages(mainP, data_DA, data_BF)

%%
speeds = [-0.5, 0, 0.5];
speed_BF = cell([1, length(speeds)]);
speed_mainP = cell([1, length(speeds)]);
for sp=1:length(speeds)
    clearvars -except mainP speed_BF speed_mainP speeds sp
    mainP.shift.val = speeds(sp);
    main_init
    speed_BF{sp} = data_BF;
    speed_mainP{sp} = copyStruct(mainP);
end
%%
if mainP.save_plots
    figure('units','normalized','position',[.2 .3 .5 .3],'Visible','off')
else
    figure;
end
m = 1;
for sp=1:length(speeds)
    m_BF = speed_BF{sp}{m};
    s_mainP = speed_mainP{sp};
    m_P = s_mainP.P;
    thetas = s_mainP.getScanGrid(s_mainP.methods_set{m});
    radius = data_DA.Radius;
%     if mainP.interp_upsample > 0
%         radius = linspace(radius(1), radius(end), mainP.interp_upsample);
%     end
    warning('off')
%     [img, Xs, Zs] = getScanConvertedImage(m_BF, ...
%         thetas, 1e3 * radius, length(thetas), length(radius), 'spline'); % Test
    [img, Xs, Zs] = getScanConvertedImage(m_BF, ...
        thetas, 1e3 * radius, 2024, 2024, 'spline');
    warning('on')
    img = db(abs(img));
    
    shifts = s_mainP.shift.getShifts(s_mainP.P);
    shifts = interp1(1:length(shifts), shifts, ...
        linspace(1, length(shifts), length(thetas)));
    max_gain = max(img(:));
    az_lim = [NaN, -NaN];
    r_lim = [NaN, -NaN];
    for p=1:length(s_mainP.pts_range)
        for t=1:length(thetas)
            s_phantom = s_mainP.shift.shiftPositions(data_phantom, shifts(t));
            p_size = 0.1; %mm
            
            beam_az = atan(thetas(t)) * s_phantom.positions(p, 3)*1e3;
            pos_az = s_phantom.positions(p, 1)*1e3;
            pos_r = s_phantom.positions(p, 3)*1e3;
            s_az = find(abs(Xs - pos_az) < p_size, 10); % pos azimuth
%             s_az = find(abs(Xs - beam_az) < p_size, 10); % beam azimuth
            s_r = find(abs(Zs - pos_r) < p_size, 10);
            if ~isempty(s_r) && ~isempty(s_az) && abs(beam_az - pos_az) < 1
                % Create circle mask (instead of square)
%                 p_r = s_r * ones(1, length(s_az));
%                 p_az = s_az * ones(1, length(s_r));
%                 img_circle = zeros(size(img));
%                 p_circle = (p_r-mean(s_r)).^2 + (p_az'-mean(s_az)).^2 <= 10;
%                 img(s_r, s_az) = p_circle * (max_gain + 3) + ...
%                      abs(p_circle - 1).*img(s_r, s_az);
%                 img(s_r, s_az) = max_gain + 3; % square mask
                az_lim = [min(az_lim(1), pos_az - p_size), ...
                    max(az_lim(2), pos_az + p_size)];
                r_lim = [min(r_lim(1), pos_r - p_size), ...
                    max(r_lim(2), pos_r + p_size)];
            end
        end
    end
    pad = 0.5; % padding in mm
%     az_lim = [max(Xs(1), az_lim(1) - pad), min(Xs(end), az_lim(2) + pad)];
    az_lim = [-1, 1];
    r_lim = [max(Zs(1), r_lim(1) - pad), min(Zs(end), r_lim(2) + pad)];
    minXs = find(Xs >= az_lim(1), 1);
    maxXs = find(Xs >= az_lim(2), 1);
    pXs = Xs(minXs:maxXs);
    minZs = find(Zs >= r_lim(1), 1);
    maxZs = find(Zs >= r_lim(end), 1);
    pZs = Zs(minZs:maxZs);
    p_img = img(minZs:maxZs, minXs:maxXs);
    p_img(~isfinite(p_img)) = -1000; % contourf doesn't handle non-finite values
%     [xGrid, zGrid] = meshgrid(linspace(pXs(1), pXs(end), 2024), ...
%         linspace(pZs(1), pZs(end), 2024));
    [xGrid, zGrid] = meshgrid(linspace(pXs(1), pXs(end), max([500, length(pXs)])), ...
        linspace(pZs(1), pZs(end), length(pZs)));
	p_img = interp2(pXs, pZs, p_img, xGrid, zGrid, 'spline', 0);
%     subplot(2,ceil(length(mainP.methods_set)/2),m);
    subplot(1, length(speeds), sp)
%     contourf(pXs, pZs, p_img, [-1000, max_gain-10, max_gain + 3], 'ShowText','on');
    contourf(xGrid, zGrid, p_img, [max_gain-50, max_gain-10, max_gain-3, max_gain + 3], 'ShowText','on');
    colormap gray
    set(gca,'YDir','Reverse')
    xlabel('azimuth [mm]');
    ylabel('range [mm]');
    title(strcat('v=', num2str(speeds(sp),2), ' m/s'));
end
if mainP.save_plots
    saveas(gcf,  mainP.outputFileName('png'), 'png')
    saveas(gcf, mainP.outputFileName('fig'), 'fig')
else
    pause
end
close
fprintf('DAS illustation finished!\n')

%% Beampattern transmission - steered response reception
clearvars -except all_peaks
mainP = MainParameters();
mainP.num_beams = 75;
mainP.medium_range = [35, 45];
mainP.P = mainP.copyP(mainP.num_beams);
mainP.shift = Shift(ShiftType.LinearSpeed, 0, mainP.num_beams, 0, 1);
mainP.shift_per_beam = true;
mainP.save_plots = false;
mainP.speckle_load = false;

% mainP.pts_azimuth = mainP.P.Tx.SinTheta .* 40;
% mainP.pts_range = cos(mainP.P.Tx.Theta) .* 40;
center_idx = ceil(length(mainP.P.Tx.Theta)/2) + 2;
% center_idx = center_idx-1:center_idx+1;
% center_idx = [center_idx, ceil(2*center_idx / 3)];
mainP.pts_azimuth = mainP.P.Tx.SinTheta(center_idx) .* 40;
mainP.pts_range = cos(mainP.P.Tx.Theta(center_idx)) .* 40;

az_val = mainP.P.Tx.SinTheta(center_idx+1) - mainP.P.Tx.SinTheta(center_idx);
mainP.pts_azimuth = mainP.P.Tx.SinTheta(center_idx) + az_val/2;
z_val = cos(mainP.P.Tx.Theta(center_idx+1)) - cos(mainP.P.Tx.Theta(center_idx));
mainP.pts_range = cos(mainP.P.Tx.Theta(center_idx)) + z_val/2;
mainP = mainP.createOutputDir();

main_init
%%
%plotBFImages(mainP, data_DA, data_BF)
for m=1:length(mainP.methods_set)
    figure;
    for p=1:length(mainP.pts_range)
        plot(data_peaks{m}{p}.beam_trajectory(1,:), ...
            data_peaks{m}{p}.beam_trajectory(3,:)); hold on
    end
    hold off; pause; close
end

%% Illustration interference
c = 1; % m/s
A = 1;
omega = 2 * pi;
D = 5; % m
k1 = sqrt(omega^2 / c^2);
k2 = - k1;

t = 1:1/4:10; % s
x = 0;
% x = 0:1/10:D; % m

s1 = A * exp(1i*omega*t - k1 * x);
s2 = A * exp(1i*omega*t - k2 * (x - D));

figure;
plot(x, real(s1), 'b')
plot(x, real(s2), 'r')

%%
linestyle_list = {'-.','--','-',':'};
markers_list = {'+','x','diamond','o'};
colors_list = {'b','r','g','k','m','y'};
save_folder = '..\images\may\extra\steered\40mm\';

for s=1:length(data_peaks{1})
    figure;
    for m=1:length(data_peaks)
        s_resp = data_peaks{m}{s}.beam_trajectory;
        s_resp(3,:) = s_resp(3,:) - max(s_resp(3,:));
        s_resp(1,:) = sin(s_resp(1,:)) * mainP.pts_range(1);
        plot(s_resp(1,:), s_resp(3,:), colors_list{m}, 'LineWidth', 2, ...
            'LineStyle', linestyle_list{m}, 'Marker', markers_list{m});
        hold on;
    end
    az = mainP.pts_azimuth(1);
    line('XData', [az az], 'YData', ylim, 'LineWidth', 2, ...
        'LineStyle', linestyle_list{2}, 'Color', colors_list{m+1});
    xlabel('Azimuth [mm]'); ylabel('Steered response gain [dB]');
    xlim([az-1, az+1]);
%     ylim([-80, -50]);
    ylim([-10, 0]);
    legend(horzcat(mainP.methods_set, 'Scatterer point azimuth'), ...
       'Location', 'South');
%     save(horzcat(save_folder, 'data_between.mat'), 'mainP', 'data_peaks', '-v7.3')
    saveas(gcf, horzcat(save_folder, 'steered_aligned_norm_40mm.fig'), 'fig')
    saveas(gcf, horzcat(save_folder, 'steered_aligned_norm_40mm.png'), 'png')
    hold off; pause; close
end

%%
clear all
mainP = MainParameters();
mainP.P.Tx.focus       = [0 0 80]*1e-3;         % Initial electonic focus
mainP.P.Rx.focus       = [0 0 80]*1e-3;         % Initial electonic focus
mainP.P.Tx.FocRad   = 80e-3;    % Focal radius
mainP.P.MinRadImage = 70e-3;     % Minimum radius for which data is recorded [m]
mainP.medium_range = [70, 90]; % mm

mainP.num_beams = 21;
mainP.pts_range = [80];
% dist_b = 0.3 * mainP.pts_range ./ floor(mainP.num_beams/2);
% mainP.pts_azimuth = dist_b .* 0.5;
mainP.pts_azimuth = [0];
mainP.NoMLA = 1;
mainP.shift = Shift(ShiftType.RadialVar, 1/8, 17, 0, 1); % Ref Shift.m
mainP.shift_per_beam = false;
% mainP.methods_set = {'DAS', 'IAA-MBMB', 'IAA-MBMB-2', 'IAA-MBMB-4'};
mainP.save_plots = true;
mainP.speckle_load = false;

if mainP.shift.type == ShiftType.RadialVar || ...
        mainP.shift.type == ShiftType.RadialCst
    % This allows to set mainP.pts_range above as radius instead.
    % This step transforms radiuses to ranges.
    mainP.pts_range = mainP.pts_range.*...
        cos(sin(mainP.pts_azimuth./mainP.pts_range));
end
mainP.P = mainP.copyP(mainP.num_beams, mainP.NoMLA);
mainP = mainP.createOutputDir();

%%

for sub=1:4
    subplot(2,2,sub);
    xlim([-0.5, 1.5]); ylim([39.5, 40.5])
%     set(gca, 'Color', [0, 0, 0])
end
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 20*0.8 20*0.4]);
% saveas(gcf, '../results/twopoints/plotsxlim/0_1_-0.6.png', 'png');
figure(1);

%%
load_folder = '..\results\main_2_1_3\speckle42_ext\2017-06-18_23.12.27\';
pts_gain = zeros(length(mainP.pts_range), length(mainP.methods_set), ...
    length(num_beams), true_shift.num_shifts);
for b=1:length(num_beams)
    for shpos=1:true_shift.num_shifts
        file_name = strcat(load_folder, 'datapeaks_', ...
            int2str(num_beams(b)), '_', int2str(shpos), '.mat');
        load(file_name)
        for m=1:length(mainP.methods_set)
            if isfield(data_peaks{m}{1}, 'peak')
                pts_gain(1, m, b, shpos) = data_peaks{m}{1}.peak(2);
            else
                pts_gain(1, m, b, shpos) = -1000;
            end
        end
    end
end

%% Pulse plot
P = Parameters();
pulse = conv(conv(P.Tx.Excitation, P.Tx.ImpulseResp), P.Rx.ImpulseResp);
pulse = pulse ./ max(pulse);
xvals = [1:length(pulse)] ./ P.fs * 10^6;

figure; plot(xvals, pulse, 'LineWidth', 3);
xlim([xvals(1), xvals(end)]); ylim([min(pulse), max(pulse)]);
xlabel('Time [\mus]', 'FontSize', 14);
ylabel('Normalized magnitude', 'FontSize', 14);
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 20*0.8 20*0.4]);

%%
file_suffix = {'0_-06', '0_06', '-45_-06', '-45_06', ...
    '45_-06', '45_06', '90_-06', '90_06'};
for f=1:8
    figure(f);
    set(gcf, 'Units', 'normalized', 'Position', [.2 .3 .5 .3]);
    set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 20*0.8 20*0.4]);
    for sub=1:4
        subplot(2,2,sub)
        set(gca, 'YTick', [39.5, 40, 40.5], 'YLim', [39.5, 40.5], ...
            'XLim', [-0.5, 0.5], 'XTick', [-0.5:0.25:0.5], ...
            'Color', [0, 0, 0])
    end
%     file_name = strcat('../motion/motion_', int2str(f), '.png');
    file_name = strcat('../motion/motion_', file_suffix{f});
    saveas(gcf, strcat(file_name, '.png'), 'png')
    saveas(gcf, strcat(file_name, '.fig'), 'fig')
end
% close all
%%
file_name = {'../DAS_steered1', '../DAS_steered2', '../DAS_steered3'};
for f=1:3
    figure(f);
    subplot(1,2,1); xlim([-1, 1])
%     set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 25*0.8 25*0.3])
    set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 30*0.8 30*0.25], ...
        'units','normalized','position',[.2 .3 .5 .3])
    saveas(gcf, strcat(file_name{f}, '.png'), 'png')
    saveas(gcf, strcat(file_name{f}, '.fig'), 'fig')
end
    
%%
file_name = {'../DAS_frame1', '../DAS_frame2', '../DAS_frame3'};
for f=1:3
    figure(f);
    set(gca, 'FontSize', 18, 'Xlim', [-0.5, 0.5], 'XTick', [-0.5:0.25:0.5], ...
        'Ylim', [39.5, 40.5], 'YTick', [39.5:0.25:40.5])
    title(''); xlabel('azimuth [mm]', 'FontSize', 26);
    ylabel('range [mm]', 'FontSize', 26);
%     set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 25*0.8 25*0.3])
    saveas(gcf, strcat(file_name{f}, '.png'), 'png')
    saveas(gcf, strcat(file_name{f}, '.fig'), 'fig')
end
   
%%
clear all
mainP = MainParameters();
% mainP.methods_set = {'DAS', 'IAA-MBMB', 'IAA-MBMB-2', 'IAA-MBMB-4'};
mainP.medium_range = [35, 45];
mainP.speckle_load = true;

rx_up = 1;
mainP.num_beams = 11*rx_up;
mainP.P = mainP.copyP(mainP.num_beams);
mainP.shift_per_beam = true;
mainP.shift = Shift(ShiftType.LinearSpeed, 0,  mainP.num_beams, 0, 1);
% mainP.shift_per_beam = false;
% mainP.shift = Shift(ShiftType.LinearSpeed, 0,  1, 0, 1);

origin_folder = mainP.save_folder;
if ~exist('speckle', 'var')
    speckle = {''}; %, '..\data\2_1_speckle_2_10-6.mat', '..\data\2_1_speckle_42_10-6.mat'};
end
mainP.pts_range = [40];
distb = mainP.pts_range(1) .* (0.3 / floor(mainP.num_beams/rx_up/2));
mainP.pts_azimuth = [0];

mainP.save_plots = true;
mainP = mainP.createOutputDir();