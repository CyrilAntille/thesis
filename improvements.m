%% Improvements - Array size
mainP = MainParameters();
mainP.num_beams = 65;
mainP.speckle_load = false;
mainP.shift_per_beam = true;
mainP.methods_set = {'DAS', 'MV', 'IAA-MBSB', 'IAA-MBMB'};
mainP.pts_range = [40, 40];
distb = mainP.pts_range(1) .* (0.3 / floor(mainP.num_beams/2));
mainP.pts_azimuth = [0, 2.5*distb];
mainP.shift = Shift(ShiftType.LinearSpeed, 0, mainP.num_beams, 0, 1);
mainP.P = mainP.copyP(mainP.num_beams, mainP.NoMLA);

mainP.save_plots = true;
mainP = mainP.createOutputDir();
origin_folder = mainP.save_folder;

%% Improvements - Array size
mainP.P.Tx.no_elements = 24;
mainP.P.Tx.ElPos = [[-(mainP.P.Tx.no_elements-1)/2:(mainP.P.Tx.no_elements-1)/2]...
    *mainP.P.Tx.pitch; zeros(2,mainP.P.Tx.no_elements)];
mainP.save_folder = strcat(origin_folder, 'array_size/');
mainP = mainP.createOutputDir();
main_2_2
clearvars -except mainP origin_folder

mainP.P.Tx.no_elements = 192;
mainP.P.Tx.ElPos = [[-(mainP.P.Tx.no_elements-1)/2:(mainP.P.Tx.no_elements-1)/2]...
    *mainP.P.Tx.pitch; zeros(2,mainP.P.Tx.no_elements)];
mainP.save_folder = strcat(origin_folder, 'array_size/');
mainP = mainP.createOutputDir();
main_2_2
clearvars -except mainP origin_folder

mainP.P.Tx.no_elements = 96;
mainP.P.Tx.ElPos = [[-(mainP.P.Tx.no_elements-1)/2:(mainP.P.Tx.no_elements-1)/2]...
    *mainP.P.Tx.pitch; zeros(2,mainP.P.Tx.no_elements)];

%% Improvements - Time-based upsampling
mainP.NoMLA = 3;
mainP.P = mainP.copyP(mainP.num_beams, mainP.NoMLA);
mainP.save_folder = strcat(origin_folder, 'time/');
mainP = mainP.createOutputDir();
main_2_2
clearvars -except mainP origin_folder

mainP.NoMLA = 5;
mainP.P = mainP.copyP(mainP.num_beams, mainP.NoMLA);
mainP.save_folder = strcat(origin_folder, 'time/');
mainP = mainP.createOutputDir();
main_2_2
clearvars -except mainP origin_folder

mainP.NoMLA = 9;
mainP.P = mainP.copyP(mainP.num_beams, mainP.NoMLA);
mainP.save_folder = strcat(origin_folder, 'time/');
mainP = mainP.createOutputDir();
main_2_2
clearvars -except mainP origin_folder

mainP.NoMLA = 1;
mainP.P = mainP.copyP(mainP.num_beams, mainP.NoMLA);

%% Improvements - Phase-based upsampling
mainP.methods_set = {'DAS', 'IAA-MBMB', 'IAA-MBMB-2', 'IAA-MBMB-4'};
mainP.save_folder = strcat(origin_folder, 'phase/');
mainP = mainP.createOutputDir();
main_2_2
clearvars -except mainP origin_folder
mainP.methods_set = {'DAS', 'MV', 'IAA-MBSB', 'IAA-MBMB'};
