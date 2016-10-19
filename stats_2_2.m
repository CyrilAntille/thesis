%% 2.2 - Motion between beams - Stats
clear all
mainP = MainParameters();
mainP.pts_range = [40, 50]; % Add a range (in mm) for each point
mainP.shift = Shift(ShiftType.RadialVar, 1/2, -1, -180, 5);
% mainP.shift = Shift(ShiftType.LinearSpeed, 1, -1, 0, 5);
mainP.num_beams = 505; % can be a single value or list of values
mainP.shift_per_beam = true;

speeds = 0:1/4:2; % Unit depends on ShiftType
mainP.save_plots = true;
mainP = mainP.createOutputDir();

%%
pts_3dB_width = zeros([length(mainP.pts_range), length(mainP.num_beams),...
    length(mainP.methods_set), length(speeds)]);
for sp=1:length(speeds)
    fprintf('Stats_2_2: Running main_2_2 with speed value: %0.2f.\n', speeds(sp));
    mainP.shift.val = speeds(sp);
    main_2_2
    for p=1:size(pts_3dB,1)
        for m=1:size(pts_3dB,2)
            for b=1:size(pts_3dB,3)
                pts_3dB_width(p,m,b,sp) = pts_3dB{p,m,b}(2);
            end
        end
    end
    clearvars -except mainP pts_3dB_width speeds
end
output_file = mainP.outputFileName(mainP.speckle_load);
fprintf('\nNSaving speed  data  into: %s\n', output_file)
save(output_file, 'mainP', 'pts_3dB_width', 'speeds', '-v7.3')

%% Plots
linestyle_list = {':','-','--','-.','-'};
markers_list = {'+','x','diamond','o','*'};

if mainP.save_plots
    figure('units','normalized','position',[.2 .3 .5 .3],'Visible','off')
else
    figure;
end

for p=1:size(pts_3dB_width, 1)
    for b=1:size(pts_3dB_width, 2)
        p1 = plot(speeds, squeeze(pts_3dB_width(p,b,:,:))', 'LineWidth', 2);
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
            im_name = strcat('speeds_p', int2str(mainP.pts_range(p)), '_',...
                int2str(mainP.num_beams(b)), '_', char(mainP.shift.type),...
                '_', int2str(mainP.shift.direction));
    %         saveas(gcf, strcat('../images/fig/', im_name, '.fig'), 'fig')
            saveas(gcf, strcat(mainP.save_folder, 'png/', im_name, '.png'), 'png')
        else
            pause
        end
    end
end
close
fprintf('Stats_2_2 finished!')
