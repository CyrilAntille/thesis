%% 2.1: Motion between frames - max scalloping loss vs range
clear all
mainP = MainParameters();
mainP.num_beams = 61;
mainP.shift = Shift(ShiftType.RadialVar, 1/2, 2, 0, 1); % Ref Shift.m
mainP.shift_per_beam = false;
mainP.methods_set = {'DAS','MV','IAA-MBSB','IAA-MBMB'};
mainP.save_plots = false;
mainP.speckle_load = true;
mainP.save_all_data = false;
mainP.normalize_bfim = false;
mainP.interp_upsample = 0;

if mainP.shift.type == ShiftType.RadialVar || ...
        mainP.shift.type == ShiftType.RadialCst
    % This allows to set mainP.pts_range above as radius instead.
    % This step transforms radiuses to ranges.
    mainP.pts_range = mainP.pts_range.*...
        cos(sin(mainP.pts_azimuth./mainP.pts_range));
end
mainP.P = mainP.copyP(mainP.num_beams);
mainP = mainP.createOutputDir();

%% Max scalloping loss
% pts_range = 35:5:55;
pts_range = 40;
mainP.pts_azimuth = [0];

max_loss = zeros(length(pts_range), length(mainP.methods_set));
for p=1:length(pts_range)
    mainP.pts_range = [pts_range(p)];
    fprintf('Main_2_1_2: Point range: %d.\n', mainP.pts_range(1));
    mainP.files_prefix = strcat(int2str(mainP.pts_range(1)), 'mm_');
    main_init
    for m=1:length(mainP.methods_set)
        min_gain = Inf; max_gain = -Inf;
        for s=1:mainP.shift.num_shifts
            sgain = data_peaks{m}{s}{1}.peak(2);
            min_gain = min([min_gain sgain]);
            max_gain = max([max_gain sgain]);
        end
        max_loss(p, m) = max_gain - min_gain;
    end
    plotBFImages(mainP, data_DA, data_BF)
    clearvars -except mainP pts_range max_loss
end
mainP.files_prefix = '';

%% Plot
linestyle_list = {'-.','--','-',':'};
markers_list = {'+','x','diamond','o'};
colors_list = {'b','r','g','k','m','c'};
if mainP.save_plots
    figure('units','normalized','position',[.2 .3 .5 .3],'Visible','off')
else
    figure;
end

pl = plot(pts_range, max_loss, 'LineWidth', 2);
for pidx=1:length(pl)
    pl(pidx).Marker = markers_list{pidx};
    pl(pidx).LineStyle = linestyle_list{pidx};
    pl(pidx).Color = colors_list{pidx};
end
grid on;
legend(mainP.methods_set, 'Location', 'best');
ylabel('Max scalloping loss [dB]');
xlabel('Scatterer point radius [mm]');
if mainP.save_plots
    prefix = mainP.files_prefix;
    mainP.files_prefix = strcat('scallop_vs_range_', mainP.files_prefix);
    saveas(gcf, mainP.outputFileName('png'), 'png')
    saveas(gcf, mainP.outputFileName('fig'), 'fig')
    mainP.files_prefix = prefix;
else
    pause;
end
close;
fprintf('Main_2_1_2 finished!')