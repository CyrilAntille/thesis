%% 2.1: Motion between frames
clear all
mainP = MainParameters();
mainP.pts_range = [40, 50];
mainP.pts_azimuth = [0, 0];
mainP.num_beams = 61;
mainP.shift = Shift(ShiftType.RadialVar, 1/2, 2, 0, 1); % Ref Shift.m
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

%% Max loss vs beams
num_beams = 61:10:211;
pts_gain = zeros(length(mainP.pts_range), length(mainP.methods_set), ...
    length(num_beams), mainP.shift.num_shifts);
for b=1:length(num_beams)
    mainP.num_beams = num_beams(b);
    fprintf('Main_2_1_3: Beams number: %d.\n', mainP.num_beams);
    mainP.P = mainP.copyP(mainP.num_beams);
    main_init
    for s=1:mainP.shift.num_shifts
        for m=1:length(mainP.methods_set)
            for p=1:length(mainP.pts_range)
                pts_gain(p, m, b, s) = data_peaks{s}{m}{p}.peak(2);
            end
        end
    end
%     plotBFImages(mainP, data_DA, data_BF)
    clearvars -except mainP num_beams pts_gain
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
    for m=1:length(mainP.methods_set)
        scallop_loss = zeros(1, length(num_beams));
        for b=1:length(num_beams)
            scallop_loss(:,b) = max(pts_gain(p, m, b, :)) - ...
                min(pts_gain(p, m, b, :));
        end
        plot(num_beams, scallop_loss, 'LineWidth', 2, 'LineStyle', ...
            linestyle_list{mod(m,length(linestyle_list))+1}, ...
            'Color', colors_list{mod(m,length(colors_list))+1}, ...
            'Marker', markers_list{mod(m,length(markers_list))+1})
        hold on;
    end

    line('XData', [num_beams(1) num_beams(end)], 'YData', [1 1], ...
        'LineWidth', 2, 'LineStyle', linestyle_list{mod(length(mainP.methods_set)+1,length(linestyle_list))+1}, ...
        'Color', colors_list{mod(length(mainP.methods_set)+1,length(colors_list))+1});
    hold off;
    legend([mainP.methods_set '1dB threshold'], 'Location', 'best');
    ylabel('Max scalloping loss [dB]');
    xlabel('Number of transmitted beams');
    if length(num_beams) > 1
        xlim([num_beams(1) num_beams(end)]) 
    end
%     ylim([0 5])
    if mainP.save_plots
        mainP.files_prefix = strcat('loss_beams_p', ...
            int2str(mainP.pts_range(p)), '_');
        output_file = mainP.outputFileName(true);
        saveas(gcf, output_file, 'png')
    else
        pause
    end
end
mainP.files_prefix = '';
close
fprintf('Main_2_1_3 finished!')