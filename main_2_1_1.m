%% 2.1: Motion between frames - loss vs shift
clear all
mainP = MainParameters();
mainP.pts_range = [55];
mainP.pts_azimuth = [0];
mainP.num_beams = 31;
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

%% Loss vs shift
pts_gain = zeros(length(mainP.pts_range), length(mainP.methods_set), ...
    mainP.shift.num_shifts);
main_init
for s=1:mainP.shift.num_shifts
    for m=1:length(mainP.methods_set)
        for p=1:length(mainP.pts_range)
            pts_gain(p, m, s) = data_peaks{m}{s}{p}.peak(2);
        end
    end
end
% plotBFImages(mainP, data_DA, data_BF)
%     clearvars -except mainP num_beams pts_gain
%%
for m=1:length(mainP.methods_set)
    maxbf = max(data_BF{1}{1}(:));
    for s=1:mainP.shift.num_shifts
        data_BF{m}{s} = data_BF{m}{s} ./ maxbf;
    end
end

data_peaks = cell([1, length(mainP.methods_set)]);
for m=1:length(mainP.methods_set)
    data_peaks{m} = cell([1, mainP.shift.num_shifts]);
    for s=1:mainP.shift.num_shifts
        data_peaks {m}{s} = computePeaksInfo(mainP, data_phantom{s}, ...
            data_DA{s}.Radius, data_BF{m}{s}, mainP.methods_set{m});
    end
end

for m=1:length(mainP.methods_set)
    for s=1:mainP.shift.num_shifts
        for p=1:length(mainP.pts_range)
            max_gain = -Inf;
            for s2=1:mainP.shift.num_shifts
                max_gain = max([max_gain data_peaks{m}{s2}{p}.peak(2)]);
            end
            pts_gain(p, m, s) = data_peaks{m}{s}{p}.peak(2) - max_gain;
        end
    end
end

%% Plot
linestyle_list = {'-.','--','-',':'};
markers_list = {'+','x','diamond','o'};
colors_list = {'b','r','g','k','m','c'};
if mainP.save_plots
    figure('units','normalized','position',[.2 .3 .5 .3],'Visible','off')
else
    figure;
end

shifts = mainP.shift.getShifts(mainP.P); % Assumes shifts in radians
shifts = sin(shifts) ./ (mainP.P.Tx.SinTheta(2) - mainP.P.Tx.SinTheta(1));
for p=1:length(mainP.pts_range)
    pl = plot(shifts, squeeze(pts_gain(p,:,:)), 'LineWidth', 2);
    for pidx=1:length(pl)
        pl(pidx).Marker = markers_list{pidx};
        pl(pidx).LineStyle = linestyle_list{pidx};
        pl(pidx).Color = colors_list{pidx};
    end
%     ylim([-20 0])
    s = floor(shifts(1));
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
%     t = strcat('Scatterer point at ', num2str(mainP.pts_range(p),0), 'mm range, ');
%     title(t)

    if mainP.save_plots
        prefix = mainP.files_prefix;
        mainP.files_prefix = strcat('loss_shift_p', ...
            int2str(mainP.pts_range(p)), '_');
        saveas(gcf, mainP.outputFileName('png'), 'png')
        saveas(gcf, mainP.outputFileName('fig'), 'fig')
        mainP.files_prefix = prefix;
        save(strcat(mainP.save_folder, mainP.files_prefix, ...
            'results_2_1_1.mat'), 'pts_gain', '-v7.3')
    else
        pause
    end
end
close
fprintf('Main_2_1_1 finished!')