mainP = MainParameters();
% mainP.methods_set = {'DAS', 'IAA-MBMB', 'IAA-MBMB-2', 'IAA-MBMB-4'};
mainP.pts_range = [40];
mainP.pts_azimuth = [0];
mainP.medium_range = [35, 45];
mainP.shift_per_beam = true;
mainP.speckle_load = false;
mainP.speckle_file = '..\data\2_1_speckle_42_10-6.mat';
mainP.save_plots = true;
mainP = mainP.createOutputDir();

origin_folder = mainP.save_folder;
speckle = {''};
% '..\data\2_1_speckle_2_10-6.mat', '..\data\2_1_speckle_42_10-6.mat'};

%% 2.1.1 Single point moving, 65 beams
% speeds = -0.6:0.1:0.6;
% mainP.num_beams = 65;
speeds = [-0.6:0.1:0.6];
% speeds = [0:0.375:3.75] ./ 3;
mainP.num_beams = 65;
mainP.NoMLA = 1;
mainP.methods_set = {'DAS', 'MV', 'IAA-MBSB', 'IAA-MBMB'};
mainP.shift = Shift(ShiftType.LinearSpeed, 0, mainP.num_beams/3, 0, 3);
mainP.P = mainP.copyP(mainP.num_beams, mainP.NoMLA);
for spkl=1:length(speckle)
    mainP.speckle_file = speckle{spkl};
    if strcmp(mainP.speckle_file, '')
        mainP.speckle_load = false;
    else
        mainP.speckle_load = true;
    end
    mainP.save_folder = strcat(origin_folder, '2.1.1/');
    mainP = mainP.createOutputDir();
    stats_2_2
    clearvars -except mainP origin_folder speckle speeds
end
fprintf('\n-------------------------------------------- 2.1.1 finished! ----------------------------------------------------------------\n')

%% 2.1.2 Single point moving in noiseless medium, 173 beams
speeds = -0.6:0.1:0.6;
mainP.num_beams = 173;
mainP.methods_set = {'DAS', 'IAA-MBMB'};
mainP.shift = Shift(ShiftType.LinearSpeed, 0, mainP.num_beams, 0, 1);
mainP.P = mainP.copyP(mainP.num_beams);
speckle = {''};
for spkl=1:length(speckle)
    mainP.speckle_file = speckle{spkl};
    if strcmp(mainP.speckle_file, '')
        mainP.speckle_load = false;
    else
        mainP.speckle_load = true;
    end
    mainP.save_folder = strcat(origin_folder, '2.1.2/');
    mainP = mainP.createOutputDir();
    stats_2_2
    clearvars -except mainP origin_folder speckle speeds
end
fprintf('\n-------------------------------------------- 2.1.2 finished! ----------------------------------------------------------------\n')

%% 2.1.3 Single point moving in noiseless medium, 191 beams
speeds = -0.6:0.1:0.6;
mainP.num_beams = 191;
mainP.methods_set = {'DAS', 'IAA-MBSB'};
mainP.shift = Shift(ShiftType.LinearSpeed, 0, mainP.num_beams, 0, 1);
mainP.P = mainP.copyP(mainP.num_beams);
speckle = {''};
for spkl=1:length(speckle)
    mainP.speckle_file = speckle{spkl};
    if strcmp(mainP.speckle_file, '')
        mainP.speckle_load = false;
    else
        mainP.speckle_load = true;
    end
    mainP.save_folder = strcat(origin_folder, '2.1.3/');
    mainP = mainP.createOutputDir();
    stats_2_2
    clearvars -except mainP origin_folder speckle speeds
end
fprintf('\n-------------------------------------------- 2.1.3 finished! ----------------------------------------------------------------\n')

%% 2.1.2 Single point moving in noiseless medium, 981 beams
% speeds = -0.6:0.1:0.6;
speeds = 0.025:0.05:0.225;
mainP.num_beams = 981;
mainP.methods_set = {'DAS', 'MV', 'IAA-MBSB', 'IAA-MBMB'};
mainP.shift = Shift(ShiftType.LinearSpeed, 0, mainP.num_beams, 0, 1);
mainP.P = mainP.copyP(mainP.num_beams);
speckle = {''};
for spkl=1:length(speckle)
    mainP.speckle_file = speckle{spkl};
    if strcmp(mainP.speckle_file, '')
        mainP.speckle_load = false;
    else        mainP.speckle_load = true;
    end
    mainP.save_folder = strcat(origin_folder, '2.1.4/');
    mainP = mainP.createOutputDir();
    stats_2_2
    clearvars -except mainP origin_folder speckle speeds
end
fprintf('\n-------------------------------------------- 2.1.4 finished! ----------------------------------------------------------------\n')
