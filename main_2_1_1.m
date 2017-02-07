%% 2.1: Motion between frames - loss vs shift
clear all
mainP = MainParameters();
mainP.pts_range = [40, 50];
mainP.pts_azimuth = [0, 0];
mainP.num_beams = 61;
mainP.shift = Shift(ShiftType.RadialVar, 1/8, 17, 0, 1); % Ref Shift.m
mainP.shift_per_beam = false;
mainP.methods_set = {'DAS','MV','IAA-MBSB','IAA-MBMB'};
mainP.save_plots = true;
mainP.speckle_load = false;
mainP.save_all_data = false;
mainP.normalize_bfim = false;

if mainP.shift.type == ShiftType.RadialVar || ...
        mainP.shift.type == ShiftType.RadialCst
    % This allows to set mainP.pts_range above as radius instead.
    % This step transforms radiuses to ranges.
    mainP.pts_range = mainP.pts_range.*...
        cos(sin(mainP.pts_azimuth./mainP.pts_range));
end
mainP.P = mainP.copyP(mainP.num_beams);
mainP = mainP.createOutputDir();

%% Loss vs shift
pts_gain = zeros(length(mainP.pts_range), length(mainP.methods_set), ...
    mainP.shift.num_shifts);
main_init
for s=1:mainP.shift.num_shifts
    for m=1:length(mainP.methods_set)
        for p=1:length(mainP.pts_range)
            pts_gain(p, m, s) = data_peaks{s}{m}{p}.peak(2);
        end
    end
end
plotBFImages(mainP, data_DA, data_BF)
%     clearvars -except mainP num_beams pts_gain

%% Plot
linestyle_list = {'-.','--','-',':'};
markers_list = {'+','x','diamond','o'};
colors_list = {'b','r','g','k','m','c'};
if mainP.save_plots
    figure('units','normalized','position',[.2 .3 .5 .3],'Visible','off')
else
    figure;
end

shifts = (0:mainP.shift.num_shifts-1) * mainP.shift.val;
for p=1:length(mainP.pts_range)
    pl = plot(shifts, squeeze(pts_gain(p,:,:)), 'LineWidth', 2);
    for pidx=1:length(pl)
        pl(pidx).Marker = markers_list{pidx};
        pl(pidx).LineStyle = linestyle_list{pidx};
        pl(pidx).Color = colors_list{pidx};
    end
    s = 0;
    while true
        if s > ceil(max(shifts))
            break
        end
        l = line('XData', [s s], 'YData', ylim, 'LineWidth', 2, ...
            'LineStyle', linestyle_list{mod(length(pl), ...
            length(linestyle_list))+1}, ...
            'Color', colors_list{mod(length(pl), length(colors_list))+1});
        s = s + 1;
    end
    legend([mainP.methods_set, 'Transmitted beams'], 'Location', 'best');
    ylabel('Scatterer point gain [dB]');
    xlabel('Shift [ratio beams separation]');
    t = strcat('Scatterer point at ', num2str(mainP.pts_range(p),0), 'mm range, ');
%     title(t)
%     ylim([42 63])

    if mainP.save_plots
        prefix = mainP.files_prefix;
        mainP.files_prefix = strcat('loss_shift_p', ...
            int2str(mainP.pts_range(p)), '_');
        saveas(gcf, mainP.outputFileName('png'), 'png')
        saveas(gcf, mainP.outputFileName('fig'), 'fig'))
        mainP.files_prefix = prefix;
    else
        pause
    end
end
close
fprintf('Main_2_1_1 finished!')