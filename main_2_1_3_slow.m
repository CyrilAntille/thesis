%% 2.1: Motion between frames
if ~exist('mainP', 'var')
    num_beams=11:10:191;
    mainP = MainParameters();
    mainP.shift_per_beam = true;
    mainP.shift = Shift(ShiftType.RadialVar, 0, 1, 0, 1);
    mainP.medium_range = [35, 45]; % mm
    mainP.methods_set = {'DAS', 'IAA-MBMB', 'IAA-MBMB-2', 'IAA-MBMB-4'};
    mainP.save_plots = true;
    mainP.speckle_load = false;
    mainP.speckle_file = '..\data\2_1_speckle_2_10-6.mat';
    mainP.NoMLA = 1;
end
mainP.P = mainP.copyP(mainP.num_beams, mainP.NoMLA);
mainP = mainP.createOutputDir();

%%
true_shift = Shift(ShiftType.RadialVar, 1/4, 4, 0, 1);
pts_gain = zeros(length(mainP.pts_range), length(mainP.methods_set), ...
    length(num_beams), true_shift.num_shifts);
init_pos = [0, 40]; % [azimuth, range] mm
fprintf('\n-------------------------------------------- START ----------------------------------------------------------------')
for b=1:length(num_beams)
    fprintf(['\n-------------------------------------------- NUM BEAMS = ', ...
        int2str(num_beams(b)), ' ----------------------------------------------------------------'])
    mainP.num_beams = num_beams(b);
    mainP.P = mainP.copyP(mainP.num_beams);
    mainP.shift = Shift(ShiftType.RadialVar, 0, mainP.num_beams, 0, 1);
    shifts_rad = true_shift.getShifts(mainP.P); % rad
    shifts_az = init_pos(2) * sin(shifts_rad) + init_pos(1); % mm
    shifts_r = init_pos(2) * cos(sin(shifts_az./init_pos(2))); % mm
    for shpos=1:true_shift.num_shifts
        fprintf(['\nShift: ', int2str(shpos), '\n'])
        mainP.pts_azimuth = shifts_az(shpos);
        mainP.pts_range = shifts_r(shpos);
        main_init
        output_file = strcat(mainP.save_folder, 'datapeaks_', ...
            int2str(num_beams(b)), '_', int2str(shpos), '.mat');
        fprintf('\nNSaving data_peaks into: %s\n', output_file)
        save(output_file, 'mainP', 'data_peaks', '-v7.3')
        for m=1:length(mainP.methods_set)
            if isfield(data_peaks{m}{1}, 'peak')
                pts_gain(1, m, b, shpos) = data_peaks{m}{1}.peak(2);
            else
                pts_gain(1, m, b, shpos) = -1000;
            end
        end
        clearvars data_DA data_BF data_peaks
    end
end
save(strcat(mainP.save_folder, 'results_2_1_3_slow.mat'), ...
    'mainP', 'num_beams', 'pts_gain', 'true_shift', 'init_pos', '-v7.3')
fprintf('\n-------------------------------------------- FINISHED ----------------------------------------------------------------\n')

%% Plots
linestyle_list = {'-','-.','--',':'};
% markers_list = {'+','x','d','o','.','s','^','>','v','<'};
markers_list = {'s','d','^','x'};
colors_list = {'b','r','g','k','m','c'};
if mainP.save_plots
    figure('units','normalized','position',[.2 .3 .5 .3],'Visible','off')
    set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 20*0.8 20*0.4])
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
    line('XData', [num_beams(1) num_beams(end)], 'YData', [1 1], ...
        'LineWidth', 3, 'LineStyle', linestyle_list{1}, ...
        'Color', colors_list{end});
    hold on;
    pl = plot(num_beams, scallop_loss', 'LineWidth', 3, 'MarkerSize', 8);
    for pidx=1:length(pl)
        pl(pidx).Marker = markers_list{pidx};
        pl(pidx).LineStyle = linestyle_list{pidx};
        pl(pidx).Color = colors_list{pidx};
    end
    hold off;
    legend(['1dB threshold', mainP.methods_set], 'Location', 'NE', 'FontSize', 14);
    ylabel('Max scalloping loss [dB]', 'FontSize', 14);
    xlabel('Number of transmitted beams', 'FontSize', 14);
    if length(num_beams) > 1
        xlim([num_beams(1) num_beams(end)]) 
    end
%     ylim([0, 5])
    if mainP.save_plots
        mainP.files_prefix = strcat('loss_beams_p', ...
            int2str(mainP.pts_range(p)), '_');
        saveas(gcf, mainP.outputFileName('png'), 'png')
        saveas(gcf, mainP.outputFileName('fig'), 'fig')
    else
        pause
    end
end
mainP.files_prefix = '';
close
fprintf('Main_2_1_3_slow finished!')