%% 2.2: Motion between beams
% clear all
if ~exist('mainP', 'var')
    mainP = MainParameters();
    mainP.pts_range = [40 40]; % Add a range (in mm) for each point
    mainP.pts_azimuth = [0 1];
    mainP.num_beams = 101; % can be a single value or list of values
    mainP.shift = Shift(ShiftType.LinearSpeed, 1, mainP.num_beams, 0, 1);
%     mainP.shift = Shift(ShiftType.RadialVar, 7/8, mainP.num_beams, 0, 1);
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
    thetaRange = m_P.Tx.Theta;
    warning('off')
    [scanConvertedImage, Xs, Zs] = getScanConvertedImage(m_BF, ...
        thetaRange, data_DA.Radius, 2024, 2024, 'spline');
    warning('on')
    img = db(abs(scanConvertedImage));

    clf(); subplot(1,2,1)
    imagesc(Xs.*1000, Zs.*1000, img)
    caxis([-25  25]);
%     xlim([-20 20])
    ylim([min(mainP.pts_range)-5 max(mainP.pts_range)+5])
    colorbar
    colormap(gray)
    xlabel('azimuth [mm]');
    ylabel('range [mm]');
    title(['Method :', mainP.methods_set{m}, ', No Beams: ', ...
        int2str(mainP.num_beams)])

    subplot(1,2,2); hold on
    bf_img = db(abs(m_BF));
    thetaRange = rad2deg(thetaRange);
%     thetaRange = sind(thetaRange) * 40; % TEST. mm
    pl_legend = {};
    y_min = Inf; y_max = -Inf; x_min = Inf; x_max = -Inf;
    for p=1:length(mainP.pts_range)
        p_range = mainP.pts_range(p); % range not constant if LateralShift
        minr = find(data_DA.Radius >= (p_range - 2)*1e-3, 1);
        maxr = find(data_DA.Radius >= (p_range + 2)*1e-3, 1);
        p_bp = max(bf_img(minr:maxr, :), [], 1);
        if p == 1 || abs(p_range - mainP.pts_range(1)) > 0.5
            plot(thetaRange, p_bp, colors_list{p}, 'LineWidth', 2, ...
                'LineStyle', linestyle_list{p}, 'Marker', markers_list{p});
        end
        p3db = data_peaks{m}{p}.peak_3db .* 1000;
        p3db = atand(p3db./p_range);
        p_y = data_peaks{m}{p}.peak(2) - 3;
        p_title = strcat('P', int2str(p), '-', ...
            int2str(data_phantom.positions(p,3)*1000), 'mm range');
        pl_legend{end+1} = p_title;
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
    xlim([x_min-1, x_max+1])
    ylim([y_min-3 y_max])
    xlabel('angle [deg]');
    ylabel('gain [dB]');
    legend(pl_legend, 'Location', 'best');
    hold off;

    if mainP.save_plots
        saveas(gcf, mainP.outputFileName('png'), 'png')
        saveas(gcf, mainP.outputFileName('fig'), 'fig')
    else
        pause
    end
end
close


%% Points contour plots
if mainP.save_plots
    set(0, 'DefaultFigureVisible', 'off')
else
    set(0, 'DefaultFigureVisible', 'on')
end
figs = [];
for p=1:length(mainP.pts_range)
    figs(end+1) = figure(p);
end

for m=1:length(mainP.methods_set)
    m_BF = data_BF{m};
    m_P = mainP.P;
    if strcmp(mainP.methods_set(m),'IAA-MBMB-Upsampled')
        m_P = mainP.copyP(mainP.upsample_number);
    end

    thetaRange = m_P.Tx.Theta;
    if strcmp(mainP.methods_set(m),'IAA-MBMB-Upsampled')
        thetaRange = linspace(m_P.Tx.Theta(1), m_P.Tx.Theta(end), ...
            mainP.upsample_number);
    end
    warning('off')
    [scanConvertedImage, Xs, Zs] = getScanConvertedImage(m_BF, ...
        thetaRange, 1e3 * data_DA.Radius, 2024, 2024, 'spline');
    warning('on')
    img = db(abs(scanConvertedImage));

    if mainP.shift_per_beam
        b_shift = Shift(mainP.shift.type, mainP.shift.val, ...
            m_P.Tx.NTheta, mainP.shift.direction);
        beam_shift = m_P.Tx.SinTheta(2) - m_P.Tx.SinTheta(1);
    else
        b_shift = mainP.shift;
    end
    shifts = b_shift.getShifts(m_P);
    max_gain = max(img(:));
    for p=1:length(mainP.pts_range)
        for s=1:b_shift.num_shifts
%             for s=floor(b_shift.num_shifts*1/3):ceil(b_shift.num_shifts*2/3)
            s_phantom = b_shift.shiftPositions(data_phantom, shifts(s));
            p_size = 0.05; %mm
%                 s_az = find(abs(Xs - s_phantom.positions(p, 1)*1e3) < p_size, 10);
            beam_az = atan(m_P.Tx.Theta(s)) * s_phantom.positions(p, 3)*1e3;
            s_az = find(abs(Xs - beam_az) < p_size, 10);
            s_r = find(abs(Zs - s_phantom.positions(p, 3)*1e3) < p_size, 10);
            if ~isempty(s_r) && ~isempty(s_az)
                % Create circle mask (instead of square)
                p_r = s_r * ones(1, length(s_az));
                p_az = s_az * ones(1, length(s_r));
                img_circle = zeros(size(img));
                p_circle = (p_r-mean(s_r)).^2 + (p_az'-mean(s_az)).^2 <= 10;
                img(s_r, s_az) = p_circle * (max_gain + 3) + ...
                     abs(p_circle - 1).*img(s_r, s_az);
            end
        end

        p_range = mainP.pts_range(p);

        if p_range == mainP.P.Tx.FocRad * 1000
            max_az = 1.5 + cosd(mainP.shift.direction);
        else
            max_az = 3 + 2 * cosd(mainP.shift.direction);
        end
        max_r = 0.5 + abs(sind(mainP.shift.direction));

        minXs = find(Xs >= -max_az, 1);
        maxXs = find(Xs >= max_az, 1);
        pXs = Xs(minXs:maxXs);
        minZs = find(Zs >= p_range-max_r, 1);
        maxZs = find(Zs >=p_range+max_r, 1);
        pZs = Zs(minZs:maxZs);
        p_img = img(minZs:maxZs, minXs:maxXs);

        set(0, 'currentfigure', figs(p));
        subplot(2,ceil(length(mainP.methods_set)/2),m);
%         contourf(pXs, pZs, p_img, [-100, pts_3dB{p,b,m}(1)-10, ...
%             pts_3dB{p,b,m}(1)-3, max_gain + 3], 'ShowText','on');
        contourf(pXs, pZs, p_img, [-100, max_gain-10, max_gain + 3], 'ShowText','on');
        set(gca,'YDir','Reverse')
        xlabel('azimuth [mm]');
        ylabel('range [mm]');
        title(mainP.methods_set{m});
    end
end
if mainP.save_plots
    output_file = mainP.outputFileName('png');
    output_file = output_file(1:end-4); % Removes .png extension
    for p=1:length(mainP.pts_range)
        p_name = strcat(output_file, int2str(mainP.pts_azimuth(p)), ...
            int2str(mainP.pts_range(p)));
        saveas(figs(p), strcat(p_name, 'png'), 'png')
        saveas(figs(p), strcat(p_name, 'fig'), 'fig')
    end
else
    pause
end
close all
set(0, 'DefaultFigureVisible', 'on')
fprintf('Main_2_2 finished!\n')