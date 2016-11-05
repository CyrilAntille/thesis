%% 2.2: Motion between beams - 2 close points
% clear all
if ~exist('mainP', 'var')
    mainP = MainParameters();
    mainP.pts_range = [40 40]; 
    mainP.pts_azimuth = [0 1];
    mainP.num_beams = 101;
    mainP.shift = Shift(ShiftType.LinearSpeed, 0.5, mainP.num_beams, 0, 1);
%     mainP.shift = Shift(ShiftType.RadialVar, 1/2, mainP.num_beams, 0, 1);
    mainP.shift_per_beam = true;
    mainP.save_plots = false;
    
    if true && (mainP.shift.type == ShiftType.RadialVar || ...
            mainP.shift.type == ShiftType.RadialCst)
        % This allows to set mainP.pts_range above as radius instead.
        % This step transforms radiuses to ranges.
        mainP.pts_range = mainP.pts_range.*...
            cos(sin(mainP.pts_azimuth./mainP.pts_range));
    end
    mainP.P = mainP.copyP(mainP.num_beams);
    mainP = mainP.createOutputDir();
end

%%
main_init
data_peaks = computePeaksInfo(mainP, data_phantom, data_DA, data_BF);

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
    m_P = mainP.P;
    if strcmp(mainP.methods_set(m),'IAA-MBMB-Upsampled')
        m_P = mainP.copyP(mainP.upsample_number);
    end
    m_pts = data_peaks{m};
    % ---- Peaks plot --------
    clf(); subplot(1,2,2); hold on
    
    legends = {};
    m_peaks = [];
    for p=1:length(m_pts)
        if isfield(m_pts{p}, 'peak')
            m_peaks = horzcat(m_peaks, m_pts{p}.peak);
        end
        if p == 1 || mainP.shift.val ~= 0 || mainP.shift.direction ~= 0 ...
                || mainP.shift.type ~= ShiftType.LinearSpeed
            plot(1e3*m_pts{p}.beam_trajectory(1,:), m_pts{p}.beam_trajectory(3,:), ...
                'LineWidth', 2, 'LineStyle', linestyle_list{p}, ...
                'Color', colors_list{p})
            legends{end+1} = strcat('P', int2str(p), '-gain');
        end
    end
    
    % If two peaks, computes dip
    if length(m_peaks) == 2 && mainP.pts_azimuth(1) ~= mainP.pts_azimuth(2)
        [pts_az_sorted, pidx] = sort(m_peaks(1,:));
        pts_az_sorted = vertcat(pts_az_sorted, m_peaks(2,pidx));
        [pts_gain_sorted, pidx] = sort(m_peaks(2,:), 'descend');
        pts_gain_sorted = vertcat(m_peaks(1,pidx), pts_gain_sorted);

        dip = [];
        for p=1:size(pts_gain_sorted, 2)
            trj = m_pts{p}.beam_trajectory;
            offset = 0.0; % mm. To avoid 'max dip' in other peak
            start_az = find(trj(1,:) >= pts_az_sorted(1,1)+offset, 1);
            stop_az = find(trj(1,:) >= pts_az_sorted(1,2)-offset, 1);
            [gain, idx] = min(trj(3,start_az:stop_az));
            if isempty(dip) || dip(2) < gain
                dip = [trj(1,start_az + idx -1), gain];
            end
        end
        if ~isempty(dip)
            line('XData', 1e3*[dip(1) dip(1)], 'YData', [dip(2) ...
                pts_gain_sorted(2,1)], 'LineWidth', 2, 'LineStyle', ...
                linestyle_list{mod(length(m_pts), ...
                length(linestyle_list))+1}, 'Color', colors_list{...
                mod(length(m_pts), length(colors_list))+1})
            legends{end+1} = strcat('Dip depth:  ', ...
                num2str(pts_gain_sorted(2,1) - dip(2), 3), ' dB');
            ylim([dip(2)-2 pts_gain_sorted(2,1)+1])
        end

        scallop_az = pts_gain_sorted(1,2);
        line('XData', 1e3*[scallop_az scallop_az], 'YData', ...
            [pts_gain_sorted(2,2) pts_gain_sorted(2,1)], 'LineWidth', 2,...
            'LineStyle', linestyle_list{mod(length(m_pts)+1, ...
            length(linestyle_list))+1}, 'Color', colors_list{...
            mod(length(m_pts)+1, length(colors_list))+1})
        legends{end+1} = strcat('Scalloping loss: ', ...
            num2str(pts_gain_sorted(2,1)-pts_gain_sorted(2,2), 3), ' dB');

        line('XData', 1e3*[pts_az_sorted(1,1), pts_az_sorted(1,2)], 'YData',...
            [pts_gain_sorted(2,1) pts_gain_sorted(2,1)],  'LineWidth', ...
            2, 'LineStyle', linestyle_list{mod(length(m_pts)+2, ...
            length(linestyle_list))+1}, 'Color', colors_list{...
            mod(length(m_pts)+2, length(colors_list))+1})
        legends{end+1} = strcat('Distance peaks:  ', ...
            num2str(1e3*(pts_az_sorted(1,2) - pts_az_sorted(1,1)), 3),' mm');
    end
    legend(legends, 'Location', 'southoutside');
    xlim([1e3*pts_az_sorted(1, 1)-2 1e3*pts_az_sorted(1,end)+2])
    xlabel('azimuth [mm]');
    ylabel('gain [dB]');
    hold off;

    % ---- BF Contour plot --------
    subplot(1,2,1);
    warning('off')
    [scanConvertedImage, Xs, Zs] = getScanConvertedImage(m_BF, ...
        m_P.Tx.Theta, data_DA.Radius, 2024, 2024, 'spline');
    warning('on')
    img = db(abs(scanConvertedImage));
    % Contour plot azimuth bounds
    az_offset = (1.0 + 1.5 * mainP.shift.val) * 1e-3; % m
    minXs = find(Xs >= pts_az_sorted(1,1) - az_offset, 1);
    maxXs = find(Xs >= pts_az_sorted(1,2) + az_offset, 1);
    if isempty(maxXs)
        maxXs = length(Xs);
    end
    pXs = Xs(minXs:maxXs);
    
    % Contour plot range bounds
    min_r = Inf; max_r = -Inf;
    for p=1:length(m_pts)
        ptrj = m_pts{p}.p_trajectory;
        pmin = min(ptrj(2,:));
        if pmin < min_r
            min_r = pmin;
        end
        pmax = max(ptrj(2,:));
        if pmax > max_r
            max_r = pmax;
        end
    end
    r_offset = 0.5 * 1e-3; % m
    minZs = find(Zs >= min_r - r_offset, 1);
    maxZs = find(Zs >= max_r + r_offset, 1);
    if isempty(maxZs)
        maxZs = length(Zs);
    end
    pZs = Zs(minZs:maxZs);
    p_img = img(minZs:maxZs, minXs:maxXs);
    
    contour_levels = [m_peaks(2,1)-10 m_peaks(2,1)-2];
    if size(m_peaks, 1) >= 2
        contour_levels = [dip(2), (pts_gain_sorted(2,1) + dip(2))/2, ...
            pts_gain_sorted(2,2)];
    end
    [c, h] = contourf(1e3*pXs, 1e3*pZs, p_img, contour_levels, ...
        'LineColor', 'k');
    caxis([dip(2)-5 max(p_img(:))])
    contourcbar
    clabel(c,'FontSize',10,'Color','r','Rotation',0,'FontWeight','bold')
    
%     imagesc(1e3*pXs, 1e3*pZs, p_img)
%     caxis([-25 25])
%     colorbar
%     colormap(gray)
    
    set(gca,'YDir','Reverse')
    xlabel('azimuth [mm]');
    ylabel('range [mm]');
    title(mainP.methods_set{m});

    if mainP.save_plots
        prefix = mainP.files_prefix;
        mainP.files_prefix = strcat(mainP.files_prefix, ...
            mainP.methods_set{m}, '_');
        output_file = mainP.outputFileName(true);
        saveas(gcf, output_file, 'png')
        mainP.files_prefix = prefix;
    else
        pause
    end
end
close

