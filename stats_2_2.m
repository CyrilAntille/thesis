%% 2.2 - Motion between beams - Stats
save_all_data = false;
enable_plots = false;
load_speckle = true;

num_beams = 101;
speeds = [0, 0.5, 1, 1.5, 2];
mainlobes_3dBwidth = cell(size(speeds));
for sp=1:length(speeds)
    if load_speckle
%         load 2_1_speckle_2_10-6.mat
        load ..\data\2_1_speckle_2_10-6.mat
    end
    shift = Shift(ShiftType.LateralSpeed, speeds(sp), -1);
    main_2_2

    mainlobes_3dBwidth{sp} = points_3dBwidth;
    clearvars -except save_all_data enable_plots ...
        load_speckle num_beams speeds mainlobes_3dBwidth
end


%% Statistics on data
for b=1:length(num_beams)
    fprintf('\nNum Beams: %d.\n', num_beams(b))
    num_m = size(mainlobes_3dBwidth{1},1);
    p1_raw = zeros(length(speeds), num_m);
    p2_raw = zeros(length(speeds), num_m);
    for s=1:length(speeds)
        for m=1:num_m
            p1_raw(s,m) = mainlobes_3dBwidth{s}{m,b}(1);
            p2_raw(s,m) = mainlobes_3dBwidth{s}{m,b}(2);
        end
    end
    figure;
    subplot(2,1,1)
    plot(speeds, p1_raw)
    p1 = plot(speeds, p1_raw, 'LineWidth', 2);
    linestyle_list = {'-.','--','-',':'};
    markers_list = {'+','x','diamond','o'};
    for pidx=1:length(p1)
        p1(pidx).Marker = markers_list{pidx};
        p1(pidx).LineStyle = linestyle_list{pidx};
    end
    legend('DAS','MV', 'IAA-MBMB','IAA-MBMB-Upsampled')
    title('Beampattern mainlobes 3dB width')
    xlabel('Speed [m/s]')
    ylabel('Mainlobe 3dB width [degrees]')

    b_rel = p1_raw;
    for m=1:num_m
        b_rel(:,m) = b_rel(:,m) ./ b_rel(1,m) * 100 - 100;
    end
    subplot(2,1,2)
    p2 = plot(speeds, b_rel, 'LineWidth', 2);
    for pidx=1:length(p2)
        p2(pidx).Marker = markers_list{pidx};
        p2(pidx).LineStyle = linestyle_list{pidx};
    end
    legend('DAS','MV', 'IAA-MBMB','IAA-MBMB-Upsampled')
    title('Mainlobes width relative increase')
    xlabel('Speed [m/s]')
    ylabel('Mainlobe width increase [%]')
    pause
end



