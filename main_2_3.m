%% 2.2: Motion between beams - 2 close points
% clear all
if ~exist('mainP', 'var')
    mainP = MainParameters();
    mainP.pts_range = [40 40]; % Add a range (in mm) for each point
    mainP.pts_azimuth = [0 4];
    mainP.num_beams = 101; % can be a single value or list of values
%     mainP.shift = Shift(ShiftType.LinearSpeed, 0.5, -1, 0, 1);
    mainP.shift = Shift(ShiftType.RadialVar, 1/2, -1, 0, 1);
    mainP.shift_per_beam = true;
    mainP.save_plots = false;
    
    mainP = mainP.createOutputDir();
end

%%
main_init

%% Plots
linestyle_list = {'-.','--','-',':'};
markers_list = {'+','x','diamond','o'};
colors_list = {'b','r','g','k','m','c'};

if mainP.save_plots
    figure('units','normalized','position',[.2 .3 .5 .3],'Visible','off')
else
    figure;
end
for m=1:length(mainP.methods_set)
    m_BF = data_BF{m};
    for b=1:length(mainP.num_beams)
        num_beams = mainP.num_beams(b);
        if strcmp(mainP.methods_set(m),'IAA-MBMB-Upsampled')
            num_beams = mainP.upsample_number;
        end
        Pb = mainP.copyP(num_beams);
        b_BF = m_BF{b};
        b_DA = data_DA{b};

        warning('off')
        [scanConvertedImage, Xs, Zs] = getScanConvertedImage(b_BF, ...
            Pb.Tx.Theta, 1e3 * b_DA.Radius, 2024, 2024, 'spline');
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
        
        bf_img = db(b_BF);
        ratio = mainP.num_beams(b) / Pb.Tx.NTheta;
        b_shift = Shift(mainP.shift.type, ratio * mainP.shift.val, ...
            Pb.Tx.NTheta, mainP.shift.direction);
        shifts = b_shift.getShifts(Pb);
        
        legends = {};
        
        pts_peaks = [];
        pts_trajectory = {};
        for p=1:length(mainP.pts_range)
            p_trajectory = [];
            for s=1:b_shift.num_shifts
                s_phantom = b_shift.shiftPositions(data_phantoms, shifts(s));
                s_az = tan(Pb.Tx.Theta(s)) * s_phantom.positions(p, 3);
                s_radius = sqrt(s_phantom.positions(p, 3)^2 + s_az^2);
                s_radius_idx = find(b_DA.Radius >= s_radius,1);
                if ~isempty(s_radius_idx)
                    p_trajectory = horzcat(p_trajectory, ...
                        [s_az * 1e3; s_radius; bf_img(s_radius_idx, s)]);
                end
            end
            pts_trajectory{end+1} =  p_trajectory;
            [pdb, paz] = findpeaks(p_trajectory(3,:), ...
                p_trajectory(1,:), 'SortStr','descend');
            if ~isempty(pts_peaks) && abs(paz(1)-pts_peaks(2,1)) <= 1e-5
                pdb = pdb(2:end); paz = paz(2:end);
            end
            pts_peaks = horzcat(pts_peaks, [pdb(1); paz(1)]);
            
            if b_shift.direction == 0 && p > 1 && b_shift.type ~= ...
                    ShiftType.RadialCst && b_shift.type ~= ShiftType.RadialVar
                continue % Pts plots overlap, so plot only one
            end
            plot(p_trajectory(1,:), p_trajectory(3,:), 'LineWidth', 2, ...
                'LineStyle', linestyle_list{p}, 'Color', colors_list{p})
            legends{end+1} = strcat('P', int2str(p), '-gain');
        end
        [pts_gain_sorted, pidx] = sort(pts_peaks(1,:), 'descend');
        pts_gain_sorted = vertcat(pts_gain_sorted, pts_peaks(2,pidx));
        
        [pts_az_sorted, pidx] = sort(pts_peaks(2,:));
        pts_az_sorted = vertcat(pts_peaks(1,pidx), pts_az_sorted);
        
        dip = [];
        for p=1:length(pts_trajectory)
            trj = pts_trajectory{p};
            offset = 0.5; % mm. To avoid 'max dip' in other peak
            start_az = find(trj(1,:) >= pts_az_sorted(2,1)+offset, 1);
            stop_az = find(trj(1,:) >= pts_az_sorted(2,2)-offset, 1);
            [gain, idx] = min(trj(3,start_az:stop_az));
            if isempty(dip) || dip(2) < gain
                dip = [gain, trj(1,start_az + idx -1)];
            end
        end
        line('XData', [dip(2) dip(2)], 'YData', [dip(1) ...
            pts_gain_sorted(1,1)], 'LineWidth', 2, 'LineStyle', ...
            linestyle_list{mod(length(mainP.pts_range), ...
            length(linestyle_list))+1}, 'Color', colors_list{...
            mod(length(mainP.pts_range), length(colors_list))+1})
        legends{end+1} = strcat('Dip depth:  ', ...
            num2str(pts_gain_sorted(1,1) - dip(2), 3), ' dB');

        scallop_az = pts_gain_sorted(2,2);
        line('XData', [scallop_az scallop_az], 'YData', ...
            [pts_gain_sorted(1,2) pts_gain_sorted(1,1)], 'LineWidth', 2,...
            'LineStyle', linestyle_list{mod(length(mainP.pts_range)+1, ...
            length(linestyle_list))+1}, 'Color', colors_list{...
            mod(length(mainP.pts_range)+1, length(colors_list))+1})
        legends{end+1} = strcat('Scalloping loss: ', ...
            num2str(pts_gain_sorted(1,1)-pts_gain_sorted(1,2), 3), ' dB');

        line('XData', [pts_az_sorted(2,1), pts_az_sorted(2,2)], 'YData',...
            [pts_gain_sorted(1,1) pts_gain_sorted(1,1)],  'LineWidth', 2, ...
            'LineStyle', linestyle_list{mod(length(mainP.pts_range)+2, ...
            length(linestyle_list))+1}, 'Color', colors_list{...
            mod(length(mainP.pts_range)+2, length(colors_list))+1})
        legends{end+1} = strcat('Distance peaks:  ', ...
            num2str(pts_az_sorted(2,2) - pts_az_sorted(2,1), 3), ' mm');
        
        legend(legends, 'Location', 'southoutside');
        ylim([dip(1)-2 pts_gain_sorted(1,1)+1])
        xlim([pts_az_sorted(2, 1)-2 pts_az_sorted(2,end)+2])
        xlabel('azimuth [mm]');
        ylabel('gain [dB]');
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

