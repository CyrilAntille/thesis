%% 2.1: Motion between frames
clear all
% clearvars -except mainP
if ~exist('mainP', 'var')
    mainP = MainParameters();
    mainP.pts_range = [40];
    mainP.pts_azimuth = [0];
    mainP.num_beams = 31;
    mainP.shift = Shift(ShiftType.RadialVar, 1/2, 2, 0, 1); % Ref Shift.m
    mainP.shift_per_beam = false;
%     mainP.shift = Shift(ShiftType.RadialVar, 1/8, 8, 0, 1); % Ref Shift.m
%     mainP.methods_set = {'DAS', 'IAA-MBMB', 'IAA-MBMB-2', 'IAA-MBMB-4'};
    mainP.save_plots = true;
    mainP.speckle_load = false;
    mainP.speckle_file = '..\data\2_1_speckle_2_10-6.mat';
    mainP.NoMLA = 1;

    if mainP.shift.type == ShiftType.RadialVar || ...
            mainP.shift.type == ShiftType.RadialCst
        % This allows to set mainP.pts_range above as radius instead.
        % This step transforms radiuses to ranges.
        mainP.pts_range = mainP.pts_range.*...
            cos(sin(mainP.pts_azimuth./mainP.pts_range));
    end
end
mainP.P = mainP.copyP(mainP.num_beams, mainP.NoMLA);
mainP = mainP.createOutputDir();

%% Max loss vs beams
% num_beams = horzcat(15:4:103, 111:4:127);
num_beams = 851:100:1051;
pts_gain = zeros(length(mainP.pts_range), length(mainP.methods_set), ...
    length(num_beams), mainP.shift.num_shifts);
for b=1:length(num_beams)
    mainP.num_beams = num_beams(b);
    fprintf('Main_2_1_3: Beams number: %d\n', num_beams(b));
    mainP.P = mainP.copyP(mainP.num_beams);
    main_init
    for s=1:mainP.shift.num_shifts
        for m=1:length(mainP.methods_set)
            for p=1:length(mainP.pts_range)
                pts_gain(p, m, b, s) = data_peaks{m}{s}{p}.peak(2);
            end
        end
    end
    scallop_loss = zeros(length(mainP.methods_set), length(num_beams));
    for m=1:length(mainP.methods_set)
        scallop_loss(m,b) = max(pts_gain(p, m, b, :)) - ...
            min(pts_gain(p, m, b, :));
    end
    scallop_loss(:,b)
%     plotBFImages(mainP, data_DA, data_BF)
    clearvars -except mainP num_beams pts_gain test
end

%% Plots
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
%         save(strcat(mainP.save_folder, mainP.files_prefix, ...
%             'results_2_1_3.mat'), 'num_beams', 'pts_gain', '-v7.3')
    else
        pause
    end
end
mainP.files_prefix = '';
close

%% Plots fit
linestyle_list = {'-.','--','-',':'};
markers_list = {'+','x','diamond','o'};
colors_list = {'b','r','g','k','m','c'};
if mainP.save_plots
    figure('units','normalized','position',[.2 .3 .5 .3],'Visible','off')
else
    figure;
end

x1db = cell([1 length(mainP.pts_range)]);
for p=1:length(mainP.pts_range)
    x1db{p} = zeros([length(mainP.methods_set) 3]);
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
%         if m == 2 % outlier rejection
%             loss_fit{m} = fit(num_beams(2:end)', ...
%                 scallop_loss(m,2:end)', 'power2');
%         end
        scatter(num_beams, scallop_loss(m,:), colors_list{m}, ...
            'Marker',  markers_list{m}, 'LineWidth', 2);
        x = 1:2000000;
        c2 = feval(loss_fit{m}, x);
        c2 = find(c2 <= 1, 1);
        cm = predint(loss_fit{m}, x);
        c1 = find(cm(:,1) <= 1, 1);
        c3 = find(cm(:,2) <= 1, 1);
        
        c1(isempty(c1)) = -1;
        c2(isempty(c2)) = -1;
        c3(isempty(c3)) = -1;
        x1db{p}(m,:) = [c1, c2, c3];
    end
%     x1db{p}
    
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
            'results_2_1_3.mat'), 'num_beams', 'pts_gain', ...
            'loss_fit', 'x1db', 'mainP', '-v7.3')
    else
        pause
    end
end
mainP.files_prefix = '';
close
fprintf('Main_2_1_3 finished!')