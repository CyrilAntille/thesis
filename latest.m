mainP = MainParameters();
% mainP.methods_set = {'DAS', 'IAA-MBMB', 'IAA-MBMB-2', 'IAA-MBMB-4'};
mainP.pts_range = [40];
mainP.pts_azimuth = [0];
mainP.medium_range = [35, 45];
mainP.shift_per_beam = true;
mainP.speckle_load = false;
mainP.save_plots = false;
mainP = mainP.createOutputDir();

rx_up = 3;
mainP.num_beams = 65*rx_up;
mainP.P = mainP.copyP(mainP.num_beams);
mainP.shift = Shift(ShiftType.LinearSpeed, 0,  mainP.num_beams, 0, 1);
distb = mainP.pts_range(1) .* (0.3 / floor(mainP.num_beams/rx_up/2));

origin_folder = mainP.save_folder;
if ~exist('speckle', 'var')
    speckle = {''}; %, '..\data\2_1_speckle_2_10-6.mat', '..\data\2_1_speckle_42_10-6.mat'};
end

%% 2.1.1: One point lateral motion
mainP.pts_range = [40];
mainP.pts_azimuth = [0];
speeds = [-0.6:0.1:0.6] / rx_up;
% speeds = [-3:0.3:3] ./ rx_up;
% speeds = [-0.6, 0, 0.6] ./ rx_up;
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
    clearvars -except mainP origin_folder rx_up distb speckle speeds
end
fprintf('\n-------------------------------------------- 2.1.1 finished! ----------------------------------------------------------------\n')

%% 2.1.2 One point multi-direction
mainP.pts_range = [40];
mainP.pts_azimuth = [0];
speeds = [-0.6, 0.6] ./ rx_up;
directions = [-45, 45, 90];
for spkl=1:length(speckle)
    mainP.speckle_file = speckle{spkl};
    if strcmp(mainP.speckle_file, '')
        mainP.speckle_load = false;
    else
        mainP.speckle_load = true;
    end
    mainP.save_folder = origin_folder;
    mainP = mainP.createOutputDir();
    dir_folder = mainP.save_folder;
    for d=1:length(directions)
        mainP.shift = Shift(ShiftType.LinearSpeed, 0,  mainP.num_beams, ...
            directions(d), 1);
        mainP.save_folder = dir_folder;
        mainP = mainP.createOutputDir();
        stats_2_2
        clearvars -except mainP origin_folder rx_up distb speckle speeds directions spkl dir_folder
    end
end
mainP.shift = Shift(ShiftType.LinearSpeed, 0,  mainP.num_beams, 0, 1);
fprintf('\n-------------------------------------------- 2.1.2 finished! ----------------------------------------------------------------\n')

%% 2.2.1: Two points lateral motion - 2*distb
mainP.pts_range = [40, 40];
mainP.pts_azimuth = [0, 2*distb];
% speeds = [-0.6:0.1:0.6] / rx_up;
speeds = [-0.6, 0, 0.6] / rx_up;
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
    clearvars -except mainP origin_folder rx_up distb speckle speeds
end
fprintf('\n-------------------------------------------- 2.2.1 finished! ----------------------------------------------------------------\n')

%% 2.2.2: Two points lateral motion - 2.5*distb
mainP.pts_range = [40, 40];
mainP.pts_azimuth = [0, 2.5*distb];
% speeds = [-0.6:0.1:0.6] / rx_up;
speeds = [-0.6, 0, 0.6] / rx_up;
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
    clearvars -except mainP origin_folder rx_up distb speckle speeds
end
fprintf('\n-------------------------------------------- 2.2.2 finished! ----------------------------------------------------------------\n')

%% 2.2.3: Two points multi-direction - 2*distb
mainP.pts_range = [40, 40];
mainP.pts_azimuth = [0, 2*distb];
% speeds = [-0.6:0.1:0.6] / rx_up;
speeds = [-0.6, 0.6] / rx_up;
directions = [-45, 45, 90];
for spkl=1:length(speckle)
    mainP.speckle_file = speckle{spkl};
    if strcmp(mainP.speckle_file, '')
        mainP.speckle_load = false;
    else
        mainP.speckle_load = true;
    end
    mainP.save_folder = origin_folder;
    mainP = mainP.createOutputDir();
    dir_folder = mainP.save_folder;
    for d=1:length(directions)
        mainP.shift = Shift(ShiftType.LinearSpeed, 0,  mainP.num_beams, ...
            directions(d), 1);
        mainP.save_folder = dir_folder;
        mainP = mainP.createOutputDir();
        stats_2_2
        clearvars -except mainP origin_folder rx_up distb speckle speeds directions spkl dir_folder
    end
end
mainP.shift = Shift(ShiftType.LinearSpeed, 0,  mainP.num_beams, 0, 1);
fprintf('\n-------------------------------------------- 2.2.3 finished! ----------------------------------------------------------------\n')

%% 2.2.4 Two points multi-direction - 2.5*distb
mainP.pts_range = [40, 40];
mainP.pts_azimuth = [0, 2.5*distb];
% speeds = [-0.6:0.1:0.6] / rx_up;
speeds = [-0.6, 0.6] / rx_up;
directions = [-45, 45, 90];
for spkl=1:length(speckle)
    mainP.speckle_file = speckle{spkl};
    if strcmp(mainP.speckle_file, '')
        mainP.speckle_load = false;
    else
        mainP.speckle_load = true;
    end
    mainP.save_folder = origin_folder;
    mainP = mainP.createOutputDir();
    dir_folder = mainP.save_folder;
    for d=1:length(directions)
        mainP.shift = Shift(ShiftType.LinearSpeed, 0,  mainP.num_beams, ...
            directions(d), 1);
        mainP.save_folder = dir_folder;
        mainP = mainP.createOutputDir();
        stats_2_2
        clearvars -except mainP origin_folder rx_up distb speckle speeds directions spkl dir_folder
    end
end
mainP.shift = Shift(ShiftType.LinearSpeed, 0,  mainP.num_beams, 0, 1);
fprintf('\n-------------------------------------------- 2.2.4 finished! ----------------------------------------------------------------\n')
