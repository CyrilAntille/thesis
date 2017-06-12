mainP = MainParameters();
mainP.pts_range = [40, 40];
mainP.shift_per_beam = true;
mainP.speckle_load = false;
mainP.speckle_file = '..\data\2_1_speckle_42_10-6.mat';
mainP.save_plots = true;
mainP = mainP.createOutputDir();

origin_folder = mainP.save_folder;
num_beams = [65];

%% 2.2.1: Two points linear idle
speeds = [0];
for b=1:length(num_beams)
    mainP.num_beams = num_beams(b);
    distb = mainP.pts_range(1) .* (0.3 / floor(mainP.num_beams/2));
    mainP.pts_azimuth = [0, 2*distb];
    mainP.shift = Shift(ShiftType.LinearSpeed, 0, mainP.num_beams, 0, 1);
    mainP.P = mainP.copyP(mainP.num_beams);
    mainP.save_folder = strcat(origin_folder, '2.2.1/', ...
        int2str(num_beams(b)), '/');
    mainP = mainP.createOutputDir();
    stats_2_2
    clearvars -except mainP origin_folder num_beams speeds
end
fprintf('\n-------------------------------------------- 2.2.1 finished! ----------------------------------------------------------------\n')


%% 2.2.2: Two points linear velocity
speeds = [-0.5, 0.5];
directions = [0, -45, 45, 90];
for b=1:length(num_beams)
    mainP.num_beams = num_beams(b);
    distb = mainP.pts_range(1) .* (0.3 / floor(mainP.num_beams/2));
    mainP.pts_azimuth = [0, 2*distb];
    mainP.P = mainP.copyP(mainP.num_beams);
    for d=1:length(directions)
        mainP.shift = Shift(ShiftType.LinearSpeed, 0,  mainP.num_beams, ...
            directions(d), 1);
        mainP.save_folder = strcat(origin_folder, '2.2.2/', ...
            int2str(num_beams(b)), '/');
        mainP = mainP.createOutputDir();
        stats_2_2
        clearvars -except mainP origin_folder num_beams speeds directions b
    end
end
fprintf('\n-------------------------------------------- 2.2.2 finished! ----------------------------------------------------------------\n')

