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
load(strcat(source_folder, scenario_folder{1}, result_file));
[p,~,b,s] = size(pts_gain);
gain_size = [p,length(scenario_folder),b,s];
methods_gain = cell([1, length(methods_set)]);
for m=1:length(methods_set)
    methods_gain{m} = zeros(gain_size);
end
for f=1:length(scenario_folder)
    load(strcat(source_folder, scenario_folder{f}, result_file));
    for m=1:length(methods_set)
        methods_gain{m}(:,f,:,:) = pts_gain(:,m,:,:);
    end
end

%% main_2_1_3 plots
for m=1:length(methods_set)
    mainP = methods_mainP{m};
    pts_gain = methods_gain{m};
    
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