%% 2.2: Motion between beams - 2 close points
% clear all
if ~exist('mainP', 'var')
    mainP = MainParameters();
    mainP.pts_range = [40 40]; % Add a range (in mm) for each point
    mainP.pts_azimuth = [0 1];
    mainP.num_beams = 101; % can be a single value or list of values
    mainP.shift = Shift(ShiftType.LinearSpeed, 1.5, -1, 0, 1);
%     mainP.shift = Shift(ShiftType.RadialVar, 7/8, -1, 0, 1);
    mainP.shift_per_beam = true;
    mainP.save_plots = false;
    
    mainP = mainP.createOutputDir();
end

%%
main_init

%% Dip measurements
dips = cell([length(mainP.num_beams), length(mainP.methods_set)]);
for m=1:length(mainP.methods_set)
    m_BF = data_BF{m};
    for b=1:length(mainP.num_beams)
        Pb = mainP.copyP(mainP.num_beams(b));
        b_BF = m_BF{b};
        b_DA = data_DA{b};

        thetaRange = Pb.Tx.Theta;
        if strcmp(mainP.methods_set(m),'IAA-MBMB-Upsampled')
            thetaRange = linspace(Pb.Tx.Theta(1), Pb.Tx.Theta(end), 500);
        end
        
        thetaRange = rad2deg(thetaRange);
        bf_img = db(b_BF);
        % For this point on, Assumes 2 points at single range
        p_range = find(b_DA.Radius >= mainP.pts_range(1)*1e-3, 1);
        bf_img = bf_img(p_range,:);
        [pdb, pdeg] = findpeaks(bf_img, thetaRange, 'SortStr','descend');
        p_3dbline = ones(size(thetaRange)) .* (pdb(1) - 3);
        
        pts_theta = sort([pdeg(1) pdeg(2)]);
        pts_gain = sort([pdb(1) pdb(2)]);
        
        [p_x, ~, ~, ~] = intersections(thetaRange, bf_img, ...
            thetaRange, p_3dbline, 1);
        dip_start_idx = find(p_x >= pts_theta(1), 1);
        dip_end_idx = find(p_x >= pts_theta(2), 1) - 1;
        dip = struct; dip.lowest_gain = nan;
        if ~isempty(dip_end_idx) && dip_start_idx > 1 && dip_end_idx > 1
            dip.peak1_gain = pts_gain(1); dip.peak1_angle = pts_theta(1);
            dip.peak1_width = p_x(dip_start_idx) - p_x(dip_start_idx-1);
            dip.peak2_gain = pts_gain(2); dip.peak2_angle = pts_theta(2);
            dip.peak2_width = p_x(dip_end_idx+1) - p_x(dip_end_idx);

            dip_start = p_x(dip_start_idx); dip_end = p_x(dip_end_idx);
            dip.start_angle = dip_start; dip.end_angle = dip_end;
            dip.gain = p_3dbline(1); dip.width = dip_end - dip_start;

            [pdb_min, pidx] = findpeaks(max(bf_img)-bf_img);
            lowest_angle_idx = pidx(thetaRange(pidx) > dip_start ...
                & thetaRange(pidx) < dip_end);
            if ~isempty(lowest_angle_idx)
                dip.lowest_angle = thetaRange(lowest_angle_idx(1));
                dip.lowest_gain = bf_img(lowest_angle_idx);
                dip.depth = max(dip.peak1_gain, dip.peak2_gain) - dip.lowest_gain;
            end
        end
        dips{b,m} = dip;
    end
end

%% Plots
linestyle_list = {'-.','--','-',':'};
markers_list = {'+','x','diamond','o'};
colors_list = {'b','r','g','k','m','y'};

if mainP.save_plots
    figure('units','normalized','position',[.2 .3 .5 .3],'Visible','off')
else
    figure;
end
for m=1:length(mainP.methods_set)
    m_BF = data_BF{m};
    for b=1:length(mainP.num_beams)
        Pb = mainP.copyP(mainP.num_beams(b));
        b_BF = m_BF{b};
        b_DA = data_DA{b};

        thetaRange = Pb.Tx.Theta;
        if strcmp(mainP.methods_set(m),'IAA-MBMB-Upsampled')
            thetaRange = linspace(Pb.Tx.Theta(1), Pb.Tx.Theta(end), ...
                mainP.upsample_number);
        end
        warning('off')
        [scanConvertedImage, Xs, Zs] = getScanConvertedImage(b_BF, ...
            thetaRange, 1e3 * b_DA.Radius, 2024, 2024, 'spline');
        warning('on')
        img = db(abs(scanConvertedImage));
        
        clf(); subplot(1,2,1)
        imagesc(Xs, Zs, img)
        caxis([-25  25]);
        xlim([-20 20])
        ylim([min(mainP.pts_range)-5 max(mainP.pts_range)+5])
        colorbar
        colormap(gray)
        xlabel('azimuth [mm]');
        ylabel('range [mm]');
        title(['Method :', mainP.methods_set{m}, ', No Beams: ', ...
            int2str(mainP.num_beams(b))])

        
        subplot(1,2,2); hold on
        p_range = find(b_DA.Radius >= mainP.pts_range(1)*1e-3, 1);
        bf_img = db(b_BF);
        bf_img = bf_img(p_range,:);
        thetaRange = rad2deg(thetaRange);
        
        plot(thetaRange, bf_img, 'LineWidth', 2);
        dip = dips{b,m};
        if ~isnan(dip.lowest_gain)
            line('XData', [dip.start_angle dip.end_angle], ...
                'YData', [dip.gain dip.gain], ...
                'LineWidth', 2, 'LineStyle', '-', 'Color', [0,0,0]+0.4)
            p_l1 = strcat('Dip 3dB width:  ', ...
                num2str(dip.width, 3), ' degrees');
            line('XData', [dip.lowest_angle dip.lowest_angle], 'YData',...
                [dip.lowest_gain, dip.lowest_gain+dip.depth], ...
                'LineWidth', 2, 'LineStyle', '--', 'Color', [0,0,0]+0.4)
            p_l2 = strcat('Dip depth:  ', num2str(dip.depth, 3), ' dB');
        
            legend({'Data', p_l1, p_l2}, 'Location', 'best');
            ylim([dip.lowest_gain-2 max([dip.peak1_gain dip.peak2_gain])])
            xlim([dip.peak1_angle-dip.peak1_width, dip.peak2_angle+dip.peak2_width])
        end 
        xlabel('angle [deg]');
        ylabel('gain [dB]');
%         legend(pl_legend, 'Location', 'best');
        hold off;
        
        
        if mainP.save_plots
            im_name = strcat(int2str(mainP.num_beams(b)), '_', ...
                char(mainP.shift.type), '_', int2str(mainP.shift.direction),...
                '_', num2str(mainP.shift.val,2), '_', mainP.methods_set{m});
    %         saveas(gcf, strcat('../images/fig/', im_name, '.fig'), 'fig')
            saveas(gcf, strcat(mainP.save_folder, 'png/', im_name, '.png'), 'png')
        else
            pause
        end
    end
end
close

