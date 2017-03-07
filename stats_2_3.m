%% Stats_2_3
clear all
mainP = MainParameters();
mainP.methods_set = {'DAS','MV','IAA-MBSB','IAA-MBMB'};
mainP.pts_range = [40];
% mainP.pts_azimuth = [0, 0.9586];
mainP.pts_azimuth = [0];
mainP.num_beams = 303; % can be a single value or list of values
% mainP.shift = Shift(ShiftType.LinearSpeed, 0, mainP.num_beams, 0, 3);
mainP.shift = Shift(ShiftType.RadialVar, 3*1/2, mainP.num_beams, 0, 3);
mainP.shift_per_beam = true;
mainP.save_plots = true;
mainP.speckle_load = false;

if true && mainP.shift.type == ShiftType.RadialVar && ...
        mainP.shift.type ~= ShiftType.RadialCst
    % This allows to set mainP.pts_range above as radius instead.
    % This step transforms radiuses to ranges.
    mainP.pts_range = mainP.pts_range.*...
        cos(sin(mainP.pts_azimuth./mainP.pts_range));
end
mainP.P = mainP.copyP(mainP.num_beams);
mainP = mainP.createOutputDir();

%% Various speeds
speeds = -0.5:1/4:0.5; % Unit depends on ShiftType
for sp=1:length(speeds)
    fprintf('Stats_2_2: Running main_2_2 with speed value: %0.2f.\n', speeds(sp));
    mainP.shift.val = speeds(sp);
    main_2_3
    clearvars -except mainP speeds
end
%% Various pts distances
% dists = 0.25:1/4:1;
% for d=1:length(dists)
%     fprintf('Stats_2_2: Running main_2_2 with distance: %0.2f.\n', dists(d));
%     mainP.pts_azimuth = [0 dists(d)];
%     mainP.files_prefix = strcat('Dist_', num2str(dists(d),3), '_');
%     main_2_3
%     clearvars -except mainP dists
% end
fprintf('Stats_2_3 terminated\n')