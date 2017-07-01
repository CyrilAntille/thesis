mainP = MainParameters();
% mainP.methods_set = {'DAS', 'IAA-MBMB', 'IAA-MBMB-2', 'IAA-MBMB-4'};
mainP.pts_range = [40];
mainP.pts_azimuth = [0];
mainP.medium_range = [35, 60];
mainP.shift_per_beam = true;
mainP.speckle_load = false;
mainP.save_plots = true;
mainP = mainP.createOutputDir();

rx_up = 3;
mainP.num_beams = 65*rx_up;
mainP.P = mainP.copyP(mainP.num_beams);
mainP.shift = Shift(ShiftType.LinearSpeed, 0,  mainP.num_beams, 0, 1);
distb = mainP.pts_range(1) .* (0.3 / floor(mainP.num_beams/2)) * rx_up;

origin_folder = mainP.save_folder;
speckle = {''}; %, '../data/2_1_speckle_2_10-6.mat', '../data/2_1_speckle_42_10-6.mat'};

%% 3.1.1: One point lateral motion - 36 mm range
mainP.pts_range = [36];
speeds = [-0.6:0.1:0.6] / rx_up;
for spkl=1:length(speckle)
    mainP.speckle_file = speckle{spkl};
    if strcmp(mainP.speckle_file, '')
        mainP.speckle_load = false;
    else
        mainP.speckle_load = true;
    end
    mainP.save_folder = origin_folder;
    mainP = mainP.createOutputDir();
    stats_2_2
    clearvars -except mainP origin_folder rx_up dist_b speckle speeds
end
fprintf('\n-------------------------------------------- 3.1.1 finished! ----------------------------------------------------------------\n')

%% 3.1.2: One point lateral motion - 55 mm range
mainP.pts_range = [55];
speeds = [-0.6:0.1:0.6] / rx_up;
for spkl=1:length(speckle)
    mainP.speckle_file = speckle{spkl};
    if strcmp(mainP.speckle_file, '')
        mainP.speckle_load = false;
    else
        mainP.speckle_load = true;
    end
    mainP.save_folder = origin_folder;
    mainP = mainP.createOutputDir();
    stats_2_2
    clearvars -except mainP origin_folder rx_up dist_b speckle speeds
end
fprintf('\n-------------------------------------------- 3.1.2 finished! ----------------------------------------------------------------\n')
