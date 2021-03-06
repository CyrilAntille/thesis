%% 2.2: Motion between beams - 2 close points
% clear all
if ~exist('mainP', 'var')
    mainP = MainParameters();
    mainP.pts_range = [40, 40];
    mainP.pts_azimuth = [0, 0.9586];
    % Note: at 101 beams for 35 degrees, perfect points dist
    % = 0.96789614739295739997530666425184 mm (4x beams dist)
    % (-> 0.24197403684823934999382666606296 [mm] beams dist)
%     mainP.methods_set = {'DAS','MV','IAA-MBSB','IAA-MBMB'};
%     mainP.methods_set = {'DAS','MV-0.0104','MV-0.5','MV-1'};
    mainP.num_beams = 505;
    mainP.shift = Shift(ShiftType.LinearSpeed, 0.5, mainP.num_beams, 0, 5);
%     mainP.shift = Shift(ShiftType.RadialVar, 0, mainP.num_beams, 0, 1);
    mainP.shift_per_beam = true;
    mainP.save_plots = true;
    mainP.speckle_load = false;
    
    if true && (mainP.shift.type == ShiftType.RadialVar || ...
            mainP.shift.type == ShiftType.RadialCst)
        % This allows to set mainP.pts_range above as radius instead.
        % This step transforms radiuses to rangesdn
        mainP.pts_range = mainP.pts_range.*...
            cos(sin(mainP.pts_azimuth./mainP.pts_range));
    end
    mainP.P = mainP.copyP(mainP.num_beams);
    mainP = mainP.createOutputDir();
end

%%
main_init
% plotBFImages(mainP, data_DA, data_BF)

%% Plots
linestyle_list = {'-','-.','--',':'};
% markers_list = {'+','x','d','o','.','s','^','>','v','<'};
markers_list = {'s','d','^','x'};
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
        plot(1e3*m_pts{p}.beam_trajectory(1,:), m_pts{p}.beam_trajectory(3,:), ...
            'LineWidth', 2, 'LineStyle', linestyle_list{p}, ...
            'Color', colors_list{p}) %, 'Marker', markers_list{p})
        if isfield(m_pts{p}, 'peak')
            m_peaks = horzcat(m_peaks, m_pts{p}.peak);
            az_range = find(m_pts{p}.beam_trajectory(1,:)>=m_pts{p}.peak(1),1);
            max_ampl = max(m_pts{p}.beam_trajectory(3,az_range-1:az_range+1));
            legends{end+1} = strcat('P', int2str(p), '-gain: ', ...
                num2str(max_ampl, 5), ' dB');
        else
            legends{end+1} = strcat('P', int2str(p), '-trajectory');
        end
    end
    if size(m_peaks) == 0
        continue
    end
    pts_az_sorted = m_peaks(:,1);
    pts_gain_sorted = m_peaks(:,1);
    % If two peaks, computes dip
    if size(m_peaks,2) == 2 && mainP.pts_azimuth(1) ~= mainP.pts_azimuth(2)
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
    % Variant 1: Image interpolation on absolute image
    [scanConvertedImage, Xs, Zs] = getScanConvertedImage(abs(m_BF), ...
        m_P.Tx.Theta, data_DA.Radius, 2024, 2024, 'spline');
    warning('on')
    img = db(scanConvertedImage);
    % Variant 2: Image interpolation on complex image (remove abs() from
    % getScanConvertedImage()).
    % img = db(abs(scanConvertedImage));
    
    % Contour plot azimuth bounds
    az_offset = (1.0 + 1.5 * mainP.shift.val) * 1e-3; % m
    minXs = find(Xs >= pts_az_sorted(1,1) - az_offset, 1);
    maxXs = find(Xs >= pts_az_sorted(1,end) + az_offset, 1);
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
    cmin = max(p_img(:)) - 10;
    if size(m_peaks, 2) >= 2
        contour_levels = [dip(2), (pts_gain_sorted(2,1) + dip(2))/2, ...
            pts_gain_sorted(2,end)];
        cmin = dip(2) - 5;
    end
    [c, h] = contourf(1e3*pXs, 1e3*pZs, p_img, contour_levels, ...
        'LineColor', 'k');
    caxis([cmin max(p_img(:))])
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
        saveas(gcf, mainP.outputFileName('png'), 'png')
        saveas(gcf, mainP.outputFileName('fig'), 'fig')
        mainP.files_prefix = prefix;
    else
        pause
    end
end
close

%% Steered response plot
figure; hold on;
for m=1:length(mainP.methods_set)
    m_pts = data_peaks{m};
    p = 1;
    plot(1e3*m_pts{p}.beam_trajectory(1,:), m_pts{p}.beam_trajectory(3,:), ...
        'LineWidth', 2, 'LineStyle', linestyle_list{m}, ...
        'Color', colors_list{m})
end
legend(mainP.methods_set, 'Location', 'best');
% xlim([1e3*pts_az_sorted(1, 1)-2 1e3*pts_az_sorted(1,end)+2])
xlabel('azimuth [mm]');
ylabel('gain [dB]');
grid on;
% hold off;
if mainP.save_plots
    prefix = mainP.files_prefix;
    saveas(gcf, mainP.outputFileName('png'), 'png')
    saveas(gcf, mainP.outputFileName('fig'), 'fig')
else
    pause
end
close all
fprintf('Main_2_3 finished!\n')