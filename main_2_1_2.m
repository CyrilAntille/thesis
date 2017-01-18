%% 2.1: Motion between frames - max scalloping loss vs range
clear all
mainP = MainParameters();
mainP.pts_range = [40, 50];
mainP.pts_azimuth = [0, 0];
mainP.num_beams = 61;
mainP.shift = Shift(ShiftType.RadialVar, 1/2, 2, 0, 1); % Ref Shift.m
mainP.shift_per_beam = false;
mainP.methods_set = {'DAS','MV','IAA-MBSB','IAA-MBMB'};
mainP.save_plots = false;
mainP.speckle_load = false;
mainP.save_all_data = false;

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
pts_range = 35:5:60;
mainP.pts_azimuth = [0];
mainP.num_beams = 61;
mainP.P = mainP.copyP(mainP.num_beams);
   
max_loss = zeros(length(pts_range), length(mainP.methods_set));
for p=1:length(pts_range)
    mainP.pts_range = [pts_range(p)];
    fprintf('Main_2_1_2: Point range: %d.\n', mainP.pts_range(1));
    main_init
    for m=1:length(mainP.methods_set)
        max_loss(p, m) = data_peaks{1}{m}{1}.peak(2) ...
            - data_peaks{2}{m}{1}.peak(2);
    end
%     plotBFImages(mainP, data_DA, data_BF)
    clearvars -except mainP pts_range max_loss
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

pl = plot(pts_range, max_loss, 'LineWidth', 2);
for pidx=1:length(pl)
    pl(pidx).Marker = markers_list{pidx};
    pl(pidx).LineStyle = linestyle_list{pidx};
    pl(pidx).Color = colors_list{pidx};
end

legend(mainP.methods_set, 'Location', 'best');
ylabel('Max scalloping loss [dB]');
xlabel('Scalloping point range [mm]');
if mainP.save_plots
    prefix = mainP.files_prefix;
    mainP.files_prefix = strcat('max_scallop_p', ...
        int2str(mainP.pts_range(p)), '_');
    output_file = mainP.outputFileName(true);
    saveas(gcf, output_file, 'png')
    mainP.files_prefix = prefix;
else
    pause;
end
close;
fprintf('Main_2_1_2 finished!')