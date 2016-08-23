%% 2.2 - Motion between beams - Stats
save_all_data = false;
show_plots = false;
save_plots = false;
load_speckle = false;

num_beams = 101;
speeds = 0:1:2;
mainlobes_3dBwidth = cell(size(speeds));
mainlobes_peak = cell(size(speeds));
for sp=1:length(speeds)
    if load_speckle
%         load 2_1_speckle_2_10-6.mat
        load ../data/2_1_speckle_2_10-6.mat
    end
    shift = Shift(ShiftType.LateralSpeed, speeds(sp), -1);
    main_2_2

    mainlobes_3dBwidth{sp} = points_3dBwidth;
    mainlobes_peak{sp} = points_peaksAmpl;
    clearvars -except save_all_data show_plots save_plots data_phantoms ...
        load_speckle num_beams speeds mainlobes_3dBwidth mainlobes_peak
end


%% Statistics on data
if show_plots
    figure('units','normalized','position',[.2 .3 .5 .5])
else
    figure('units','normalized','position',[.2 .3 .5 .3],'Visible','off')
end

for b=1:length(num_beams)
    fprintf('\nNum Beams: %d.\n', num_beams(b))
    for p=1:2
        num_m = size(mainlobes_3dBwidth{1},1);

        p_raw = zeros(length(speeds), num_m);
        for s=1:length(speeds)
            for m=1:num_m
                p_raw(s,m) = mainlobes_3dBwidth{s}{m,b}(p);
            end
        end
        subplot(2,1,1)
        p1 = plot(speeds, p_raw, 'LineWidth', 2);
        linestyle_list = {':','-','--','-.','-'};
        markers_list = {'+','x','diamond','o','*'};
        for pidx=1:length(p1)
            p1(pidx).Marker = markers_list{pidx};
            p1(pidx).LineStyle = linestyle_list{pidx};
        end
        legend({'DAS','MV', 'IAA-MBMB','IAA-MBMB-Upsampled'}, 'Location', 'best')
        title('Beampattern mainlobes 3dB width')
        xlabel('Speed [m/s]')
        ylabel('Mainlobe 3dB width [degrees]')

        b_rel = p_raw;
        for m=1:num_m
            b_rel(:,m) = b_rel(:,m) ./ b_rel(1,m) * 100 - 100;
        end
        subplot(2,1,2)
        p2 = plot(speeds, b_rel, 'LineWidth', 2);
        for pidx=1:length(p2)
            p2(pidx).Marker = markers_list{pidx};
            p2(pidx).LineStyle = linestyle_list{pidx};
        end
        legend({'DAS','MV', 'IAA-MBMB','IAA-MBMB-Upsampled'}, 'Location', 'best')
        title('Mainlobes width relative increase')
        xlabel('Speed [m/s]')
        ylabel('Mainlobe width increase [%]')
        if save_plots
            p_range = round(data_phantoms.positions(p,3) * 1000);
            im_name = strcat('stats_', int2str(num_beams(b)), '_', ...
                int2str(p_range), '_', num2str(speeds(1), 2), ...
                '_', num2str(speeds(end), 2));
            saveas(gcf, strcat('../images/fig/', im_name, '.fig'), 'fig')
            saveas(gcf, strcat('../images/png/', im_name, '.png'), 'png')
        end
        if show_plots
            pause
        end
    end
end
close

%% 
if show_plots
    figure('units','normalized','position',[.2 .3 .5 .5])
else
    figure('units','normalized','position',[.2 .3 .5 .3],'Visible','off')
end

for b=1:length(num_beams)
    fprintf('\nNum Beams: %d.\n', num_beams(b))
    for p=1:2
        num_m = size(mainlobes_3dBwidth{1},1);

        p_raw = zeros(length(speeds), num_m);
        for s=1:length(speeds)
            for m=1:num_m
                p_raw(s,m) = mainlobes_peak{s}{m,b}(p);
            end
        end
        subplot(2,1,1)
        p1 = plot(speeds, p_raw, 'LineWidth', 2);
        for pidx=1:length(p1)
            p1(pidx).Marker = markers_list{pidx};
            p1(pidx).LineStyle = linestyle_list{pidx};
        end
        legend({'DAS','MV', 'IAA-MBMB','IAA-MBMB-Upsampled'}, 'Location', 'best')
        title('Peak amplitude')
        xlabel('Speed [m/s]')
        ylabel('Peak amplitude [dB]')

        b_rel = p_raw;
        for m=1:num_m
            b_rel(:,m) = b_rel(:,m) ./ b_rel(1,m) * 100 - 100;
        end
        subplot(2,1,2)
        p2 = plot(speeds, b_rel, 'LineWidth', 2);
        for pidx=1:length(p2)
            p2(pidx).Marker = markers_list{pidx};
            p2(pidx).LineStyle = linestyle_list{pidx};
        end
        legend({'DAS','MV', 'IAA-MBMB','IAA-MBMB-Upsampled'}, 'Location', 'best')
        title('Mainlobes width relative increase')
        xlabel('Speed [m/s]')
        ylabel('Peak amplitude increase [%]')
        if save_plots
            p_range = round(data_phantoms.positions(p,3) * 1000);
            im_name = strcat('maxPeaks_', int2str(num_beams(b)), '_', ...
                int2str(p_range), '_', num2str(speeds(1), 2), ...
                '_', num2str(speeds(end), 2));
            saveas(gcf, strcat('../images/fig/', im_name, '.fig'), 'fig')
            saveas(gcf, strcat('../images/png/', im_name, '.png'), 'png')
        end
        if show_plots
            pause
        end
    end
end
close

fprintf('Stats_2_2 finished!')
