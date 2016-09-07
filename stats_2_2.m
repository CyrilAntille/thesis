%% 2.2 - Motion between beams - Stats
clear all
mainP = MainParameters();
mainP.pts_theta = [0, 0]; % Add a theta (in degrees) for each point
mainP.pts_range = [40*1e-3, 50*1e-3]; % Add a range (in m) for each point
mainP.shift = Shift(ShiftType.LateralSpeed, 1, -1); % Ref Shift.m
mainP.num_beams = 101; % can be a single value or list of values
mainP.shift_per_beam = true;

speeds = 0:1:2;
all_3dB_pts = cell(size(speeds));
all_3dB_inters = cell(size(speeds));
for sp=1:length(speeds)
    mainP.shift = Shift(ShiftType.LateralSpeed, speeds(sp), -1);
    main_init

    %3dB width computations
    max_peak_DAS = 0;
    pts_3dB = cell([length(mainP.pts_range), ...
        length(mainP.methods_set), length(mainP.num_beams)]);
    inters_3dB = cell([length(mainP.pts_range), ...
        length(mainP.methods_set), length(mainP.num_beams)]);
    for m=1:length(mainP.methods_set)
        m_BF = data_BF{m};
        for b=1:length(mainP.num_beams)
            Pb = mainP.copyP(mainP.num_beams(b));
            b_BF = m_BF{b};
            b_DA = data_DA{b};
            if m == 1
                max_peak_DAS = max(db(abs(b_BF(:))));
            end

            thetaRange = Pb.Tx.Theta;
            if strcmp(mainP.methods_set(m),'IAA-MBMB-Upsampled')
                thetaRange = linspace(Pb.Tx.Theta(1), Pb.Tx.Theta(end), 500);
            end

            thetaRange = rad2deg(thetaRange);
            bf_img = db(b_BF)  - max_peak_DAS;
            for p=1:length(mainP.pts_theta)
                p_range = mainP.pts_range(p); % range not constant if LateralShift
                minr = find(b_DA.Radius >= p_range - 2*1e-3, 1);
                maxr = find(b_DA.Radius >= p_range + 2*1e-3, 1);
                p_bp = max(bf_img(minr:maxr, :), [], 1);
                [pdb, pdeg] = findpeaks(p_bp, thetaRange, ...
                    'SortStr','descend');
                p_3dbline = ones(size(thetaRange)) .* (pdb(1) - 3);
                [p_x, p_y, ~, ~] = intersections(thetaRange, p_bp, ...
                    thetaRange, p_3dbline, 1);
                if length(p_x) < 2
                    p_width = NaN;
                else
                    p_width = p_x(2) - p_x(1);
                end
                pts_3dB{p,m,b} = [pdb(1), p_width];
                inters_3dB{p,m,b} = p_x;
            end
        end
    end
    all_3dB_pts{sp} = pts_3dB;
    all_3dB_inters{sp} = inters_3dB;
    clearvars -except all_3dB_pts all_3dB_inters mainP speeds
end

%% Plots
figure;
linestyle_list = {':','-','--','-.','-'};
markers_list = {'+','x','diamond','o','*'};
for b=1:length(num_beams)
    fprintf('\nNum Beams: %d.\n', num_beams(b))
    for p=1:length(mainP.pts_theta)
        num_m = size(mainlobes_3dBwidth{1},1);

        p_raw = zeros(length(speeds), num_m);
        for s=1:length(speeds)
            for m=1:num_m
                p_raw(s,m) = mainlobes_3dBwidth{s}{m,b}(p);
            end
        end
        subplot(2,1,1)
        p1 = plot(speeds, p_raw, 'LineWidth', 2);
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
