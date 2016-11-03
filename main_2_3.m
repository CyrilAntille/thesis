%% 2.2: Motion between beams - 2 close points
% clear all
if ~exist('mainP', 'var')
    mainP = MainParameters();
    mainP.pts_range = [40 40]; 
    mainP.pts_azimuth = [0 1.1];
    mainP.num_beams = 101; % can be a single value or list of values
    mainP.shift = Shift(ShiftType.LinearSpeed, 0, -1, 0, 1);
%     mainP.shift = Shift(ShiftType.RadialVar, 1/2, -1, 0, 1);
    mainP.shift_per_beam = true;
    mainP.save_plots = false;
    
    if true && (mainP.shift.type == ShiftType.RadialVar || ...
            mainP.shift.type == ShiftType.RadialCst)
        % This allows to set mainP.pts_range above as radius instead.
        % This step transforms radiuses to ranges.
        mainP.pts_range = mainP.pts_range.*...
            cos(sin(mainP.pts_azimuth./mainP.pts_range));
    end
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
    num_beams = mainP.num_beams(b);
    if strcmp(mainP.methods_set(m),'IAA-MBMB-Upsampled')
        num_beams = mainP.upsample_number;
    end
    Pb = mainP.copyP(num_beams);
    b_BF = data_BF{m}{b};
    b_DA = data_DA{b};

    % ---- Peaks plot --------
    clf(); subplot(1,2,2); hold on

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
        pt_az = [];
        for s=1:b_shift.num_shifts
            s_phantom = b_shift.shiftPositions(data_phantoms, shifts(s));
            pt_az(end+1) = s_phantom.positions(p, 1) * 1e3; % mm

            s_az = tan(Pb.Tx.Theta(s)) * s_phantom.positions(p, 3);
            s_radius = sqrt(s_phantom.positions(p, 3)^2 + s_az^2);
            s_radius_idx = find(b_DA.Radius >= s_radius,1);
            if ~isempty(s_radius_idx)
                p_trajectory = horzcat(p_trajectory, ...
                    [s_az * 1e3; s_phantom.positions(p, 3) * 1e3; ...
                    bf_img(s_radius_idx, s)]);
            end
        end
        pts_trajectory{end+1} =  p_trajectory;
        [pdb, paz] = findpeaks(p_trajectory(3,:), ...
            p_trajectory(1,:), 'SortStr','descend');

        % Expected point azimuth when beam hits it
        diff_az = abs(pt_az - p_trajectory(1,:));
        [~, az_idx] = min(diff_az);
        pt_center = p_trajectory(1, az_idx);
        % Finds nearest peak
        [~, paz_idx] = min(abs(paz - pt_center));
        pts_peaks = horzcat(pts_peaks, [pdb(paz_idx); paz(paz_idx)]);

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
        offset = 0.0; % mm. To avoid 'max dip' in other peak
        start_az = find(trj(1,:) >= pts_az_sorted(2,1)+offset, 1);
        stop_az = find(trj(1,:) >= pts_az_sorted(2,2)-offset, 1);
        [gain, idx] = min(trj(3,start_az:stop_az));
        if isempty(dip) || dip(2) < gain
            dip = [gain, trj(1,start_az + idx -1)];
        end
    end
    if ~isempty(dip)
        line('XData', [dip(2) dip(2)], 'YData', [dip(1) ...
            pts_gain_sorted(1,1)], 'LineWidth', 2, 'LineStyle', ...
            linestyle_list{mod(length(mainP.pts_range), ...
            length(linestyle_list))+1}, 'Color', colors_list{...
            mod(length(mainP.pts_range), length(colors_list))+1})
        legends{end+1} = strcat('Dip depth:  ', ...
            num2str(pts_gain_sorted(1,1) - dip(1), 3), ' dB');
        ylim([dip(1)-2 pts_gain_sorted(1,1)+1])
    end

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
    xlim([pts_az_sorted(2, 1)-2 pts_az_sorted(2,end)+2])
    xlabel('azimuth [mm]');
    ylabel('gain [dB]');
    hold off;

    % ---- BF Contour plot --------
    subplot(1,2,1);
    warning('off')
    [scanConvertedImage, Xs, Zs] = getScanConvertedImage(b_BF, ...
        Pb.Tx.Theta, 1e3 * b_DA.Radius, 2024, 2024, 'spline');
    warning('on')
    img = db(abs(scanConvertedImage));
    % Contour plot azimuth bounds
    az_offset = 1.0; % mm
    minXs = find(Xs >= pts_az_sorted(2,1) - az_offset, 1);
    maxXs = find(Xs >= pts_az_sorted(2,2) + az_offset, 1);
    if isempty(maxXs)
        maxXs = length(Xs);
    end
    pXs = Xs(minXs:maxXs);
    
    % Contour plot range bounds
    min_r = Inf; max_r = -Inf;
    for p=1:length(pts_trajectory)
        ptrj = pts_trajectory{p};
        pmin = min(ptrj(2,:));
        if pmin < min_r
            min_r = pmin;
        end
        pmax = max(ptrj(2,:));
        if pmax > max_r
            max_r = pmax;
        end
    end
    r_offset = 0.5; % mm
    minZs = find(Zs >= min_r - r_offset, 1);
    maxZs = find(Zs >= max_r + r_offset, 1);
    if isempty(maxZs)
        maxZs = length(Zs);
    end
    pZs = Zs(minZs:maxZs);
    p_img = img(minZs:maxZs, minXs:maxXs);
    
%     [c, h] = contourf(pXs, pZs, p_img, [-100, dip(1)-1, pts_gain_sorted(1,2)-1, ...
%         pts_gain_sorted(1,1)-1], 'LineColor', 'k', 'ShowText','on');
%     caxis([dip(1)-10 max(p_img(:))])
%     contourcbar
    
    imagesc(pXs, pZs, p_img)
    caxis([-25 25])
    colorbar
    colormap(gray)
    
    set(gca,'YDir','Reverse')
    xlabel('azimuth [mm]');
    ylabel('range [mm]');
    title(mainP.methods_set{m});

    if mainP.save_plots
        im_name = strcat(mainP.files_prefix, int2str(mainP.num_beams(b)), '_', ...
            char(mainP.shift.type), '_', int2str(mainP.shift.direction),...
            '_', num2str(mainP.shift.val,2), '_', mainP.methods_set{m});
%         saveas(gcf, strcat('../images/fig/', im_name, '.fig'), 'fig')
        saveas(gcf, strcat(mainP.save_folder, 'png/', im_name, '.png'), 'png')
    else
        pause
    end
end
close

