%% 2.2 - Motion between beams - Stats
clear all
mainP = MainParameters();
mainP.pts_range = [40]; % Add a range (in mm) for each point
mainP.pts_azimuth = [0]; % Add a range (in mm) for each point
mainP.num_beams = 3*71; % can be a single value or list of values
mainP.shift = Shift(ShiftType.RadialVar, 0,  mainP.num_beams, 0, 3);
% mainP.shift = Shift(ShiftType.LinearSpeed, 0,  mainP.num_beams, 0, 3);
mainP.shift_per_beam = true;
mainP.speckle_load = false;

if (mainP.shift.type == ShiftType.RadialVar || ...
        mainP.shift.type == ShiftType.RadialCst)
    % This allows to set mainP.pts_range above as radius instead.
    % This step transforms radiuses to ranges.
    mainP.pts_range = mainP.pts_range.*...
        cos(sin(mainP.pts_azimuth./mainP.pts_range));
end

mainP.save_plots = true;
mainP.P = mainP.copyP(mainP.num_beams);
mainP = mainP.createOutputDir();

%%
speeds = -0.5:1/4:0.5; % Unit depends on ShiftType
pts_3dB_width = zeros([length(mainP.pts_range),...
    length(mainP.methods_set), length(speeds)]);
for sp=1:length(speeds)
    fprintf('Stats_2_2: Running main_2_2 with speed value: %0.2f.\n', speeds(sp));
    mainP.shift.val = speeds(sp);
    main_2_2
    for p=1:size(data_peaks,1)
        for m=1:size(data_peaks,2)
            if isfield(data_peaks{m}{p}, 'peak_3db')
                pts_3dB_width(p,m,sp) = data_peaks{m}{p}.peak_3db(2);
            end
        end
    end
    clearvars -except mainP pts_3dB_width speeds
end
mainP.files_prefix = 'speeds_';
fprintf('\nNSaving speed  data  into: %s\n', mainP.outputFileName('mat'))
save(mainP.outputFileName('mat'), 'mainP', 'pts_3dB_width', 'speeds', '-v7.3')

%% Plots
linestyle_list = {':','-','--','-.','-'};
markers_list = {'+','x','diamond','o','*'};

if mainP.save_plots
    figure('units','normalized','position',[.2 .3 .5 .3],'Visible','off')
else
    figure;
end

for p=1:size(pts_3dB_width, 1)
    p1 = plot(speeds, squeeze(pts_3dB_width(p,:,:))', 'LineWidth', 2);
    for pidx=1:length(p1)
        p1(pidx).Marker = markers_list{pidx};
        p1(pidx).LineStyle = linestyle_list{pidx};
    end
    legend(mainP.methods_set, 'Location', 'best')
    title('Beampattern mainlobes 3dB width')
    x_unit = ShiftType.getShiftTypeUnit(mainP.shift.type);
    xlabel(strcat('Speed [', x_unit, ']'))
    ylabel('Mainlobe 3dB width [degrees]')
    if mainP.save_plots
        mainP.files_prefix = strcat('speeds_p', int2str(p), '_');
        saveas(gcf, mainP.outputFileName('png'), 'png')
        saveas(gcf, mainP.outputFileName('fig'), 'fig')
        mainP.files_prefix = '';
    else
        pause
    end
end
close
fprintf('Stats_2_2 finished!')
