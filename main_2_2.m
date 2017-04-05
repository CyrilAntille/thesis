%% 2.2: Motion between beams
% clear all
if ~exist('mainP', 'var')
    mainP = MainParameters();
    mainP.pts_range = [40]; % Add a range (in mm) for each point
    mainP.pts_azimuth = [0];
    mainP.num_beams = 71;
    % Note: Azimuth distance at 40 mm range: 23.6416 mm
    % 0.2 ms per beam -> beams 'velocity' = 1.665 m/s
%     mainP.shift = Shift(ShiftType.LinearSpeed, 0.5, mainP.num_beams, 45, 3);
    mainP.shift = Shift(ShiftType.LinearSpeed, 0.5, mainP.num_beams, 0, 1);
    mainP.shift_per_beam = true;
    mainP.save_plots = false;
    mainP.speckle_load = true;
    mainP.interp_upsample = 0;
    mainP.normalize_bfim = false;
    
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
plotBFImages(mainP, data_DA, data_BF)

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
    m_P = mainP.P;
    if strcmp(mainP.methods_set(m),'IAA-MBMB-Upsampled')
        m_P = mainP.copyP(mainP.upsample_number);
    end
    thetas = m_P.Tx.Theta;
    if mainP.interp_upsample > 0
        thetas = linspace(mainP.P.Tx.Theta(1), mainP.P.Tx.Theta(end), mainP.interp_upsample);
    end
    radius = data_DA.Radius;
    if mainP.interp_upsample > 0
        radius = linspace(radius(1), radius(end), mainP.interp_upsample);
    end
    warning('off')
    [scanConvertedImage, Xs, Zs] = getScanConvertedImage(m_BF, ...
        thetas, 1e3 * radius, 2024, 2024, 'spline');
    warning('on')
    img = db(abs(scanConvertedImage));

    clf(); subplot(1,2,1)
    imagesc(Xs, Zs, img)
    if mainP.normalize_bfim && (mainP.speckle_create || mainP.speckle_load)
        caxis([-25  25]);
    end
    ylim([min(mainP.pts_range)-3 max(mainP.pts_range)+3])
    colorbar
    colormap(gray)
    xlabel('azimuth [mm]');
    ylabel('range [mm]');
    title(['Method :', mainP.methods_set{m}, ', No Beams: ', ...
        int2str(mainP.num_beams)])

    subplot(1,2,2); hold on
    img = db(abs(m_BF));
    thetaRange = rad2deg(thetas);
    pl_legend = {};
    y_min = Inf; y_max = -Inf; x_min = Inf; x_max = -Inf;
    for p=1:length(mainP.pts_range)
        p_range = mainP.pts_range(p); % range not constant if LateralShift
        minr = find(radius >= (p_range - 2)*1e-3, 1);
        maxr = find(radius >= (p_range + 2)*1e-3, 1);
        p_bp = max(img(minr:maxr, :), [], 1);
        if p == 1 || abs(p_range - mainP.pts_range(1)) > 0.5
            plot(thetas, p_bp, colors_list{p}, 'LineWidth', 2, ...
                'LineStyle', linestyle_list{p}, 'Marker', markers_list{p});
            p_title = strcat('P', int2str(p), '-', ...
                int2str(data_phantom.positions(p,3)*1000), 'mm range');
            pl_legend{end+1} = p_title;
        end
        if isfield(data_peaks{m}{p}, 'peak') && ...
                isfield(data_peaks{m}{p}, 'peak_3db')
            p_y = data_peaks{m}{p}.peak(2) - 3;
            p3db = data_peaks{m}{p}.peak_3db .* 1000;
            p3db = atand(p3db./p_range);
            if length(p3db) >= 3
                line('XData', p3db(1:2), 'YData', [p_y p_y], ...
                    'LineWidth', 2, 'LineStyle', '-', 'Color', [0,0,0]+0.4)
                p_l = strcat('P', int2str(p), '-3dB width= ', ...
                    num2str(p3db(3), 3), ' [deg]');

                pl_legend{end+1} = p_l;
                x_min = min([x_min, p3db(1)]);
                x_max = max([x_max, p3db(2)]);
                y_min = min([y_min, p_y]);
                y_max = max([y_max, p_y+3]);
            end
        end
    end
    if x_min < Inf && x_max > -Inf
        xlim([x_min-1, x_max+1])
    end
    if y_min < Inf && y_max > -Inf
        ylim([y_min-3 y_max])
    end
    xlabel('angle [deg]');
    ylabel('gain [dB]');
    legend(pl_legend, 'Location', 'best');
    hold off;

    if mainP.save_plots
        mainP.files_prefix = strcat(mainP.methods_set{m}, '_');
        saveas(gcf, mainP.outputFileName('png'), 'png')
        saveas(gcf, mainP.outputFileName('fig'), 'fig')
        mainP.files_prefix = '';
    else
        pause
    end
end
close


%% Points contour plots
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

    thetas = m_P.Tx.Theta;
    if mainP.interp_upsample > 0
        thetas = linspace(mainP.P.Tx.Theta(1), mainP.P.Tx.Theta(end), mainP.interp_upsample);
    end
    radius = data_DA.Radius;
    if mainP.interp_upsample > 0
        radius = linspace(radius(1), radius(end), mainP.interp_upsample);
    end
    warning('off')
    [img, Xs, Zs] = getScanConvertedImage(m_BF, ...
        thetas, 1e3 * radius, 2024, 2024, 'spline');
    warning('on')

    if mainP.shift_per_beam
        b_shift = Shift(mainP.shift.type, mainP.shift.val, ...
            m_P.Tx.NTheta, mainP.shift.direction);
        beam_shift = m_P.Tx.SinTheta(2) - m_P.Tx.SinTheta(1);
    else
        b_shift = mainP.shift;
    end
    shifts = b_shift.getShifts(m_P);
    max_gain = max(img(:));
    az_lim = [NaN, -NaN];
    r_lim = [NaN, -NaN];
    for p=1:length(mainP.pts_range)
        for s=1:b_shift.num_shifts
            s_phantom = b_shift.shiftPositions(data_phantom, shifts(s));
            p_size = 0.05; %mm
            
            beam_az = atan(thetas(s)) * s_phantom.positions(p, 3)*1e3;
            pos_az = s_phantom.positions(p, 1)*1e3;
            pos_r = s_phantom.positions(p, 3)*1e3;
            s_az = find(abs(Xs - pos_az) < p_size, 10); % beam azimuth
%             s_az = find(abs(Xs - beam_az) < p_size, 10); % pos azimuth
            s_r = find(abs(Zs - pos_r) < p_size, 10);
            if ~isempty(s_r) && ~isempty(s_az) && abs(beam_az - pos_az) < 1
                % Create circle mask (instead of square)
                p_r = s_r * ones(1, length(s_az));
                p_az = s_az * ones(1, length(s_r));
                img_circle = zeros(size(img));
                p_circle = (p_r-mean(s_r)).^2 + (p_az'-mean(s_az)).^2 <= 10;
                img(s_r, s_az) = p_circle * (max_gain + 3) + ...
                     abs(p_circle - 1).*img(s_r, s_az);
                az_lim = [min(az_lim(1), pos_az - p_size), ...
                    max(az_lim(2), pos_az + p_size)];
                r_lim = [min(r_lim(1), pos_r - p_size), ...
                    max(r_lim(2), pos_r + p_size)];
            end
        end
    end

    pad = 0.4; % padding in mm
    az_lim = [max(Xs(1), az_lim(1) - pad), min(Xs(end), az_lim(2) + pad)];
    r_lim = [max(Zs(1), r_lim(1) - pad), min(Zs(end), r_lim(2) + pad)];
    minXs = find(Xs >= az_lim(1), 1);
    maxXs = find(Xs >= az_lim(2), 1);
    pXs = Xs(minXs:maxXs);
    minZs = find(Zs >= r_lim(1), 1);
    maxZs = find(Zs >= r_lim(end), 1);
    pZs = Zs(minZs:maxZs);
    p_img = img(minZs:maxZs, minXs:maxXs);
    p_img(p_img < -200) = -200; % contourf doesn't handle non-finite values

    subplot(2,ceil(length(mainP.methods_set)/2),m);
    contourf(pXs, pZs, p_img, [-100, max_gain-10, max_gain + 3], 'ShowText','on');
    colormap parula
    set(gca,'YDir','Reverse')
    xlabel('azimuth [mm]');
    ylabel('range [mm]');
    title(mainP.methods_set{m});
end
if mainP.save_plots
    saveas(gcf,  mainP.outputFileName('png'), 'png')
    saveas(gcf, mainP.outputFileName('fig'), 'fig')
else
    pause
end
close
fprintf('Main_2_2 finished!\n')