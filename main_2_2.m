%% 2.2: Motion between beams
% clear all
if ~exist('mainP', 'var')
    mainP = MainParameters();
    mainP.num_beams = 11;
    mainP.NoMLA = 1;
    mainP.pts_range = [40, 40];
    distb = mainP.pts_range(1) .* (0.3 / floor(mainP.num_beams/2));
    mainP.pts_azimuth = [0, 2*distb];
    mainP.methods_set = {'DAS', 'IAA-MBMB', 'MV', 'IAA-MBSB'};
%     mainP.methods_set = {'DAS', 'IAA-MBMB', 'IAA-MBMB-2', 'IAA-MBMB-5'};
    mainP.shift = Shift(ShiftType.LinearSpeed, 0, mainP.num_beams, 0, 1);
    mainP.shift_per_beam = true;
    mainP.save_plots = false;
    mainP.speckle_load = false;
    mainP.normalize_bfim = false;
    mainP.norm_variant = 2;
    mainP.P = mainP.copyP(mainP.num_beams, mainP.NoMLA);
    mainP = mainP.createOutputDir();
end

%%
main_init
% plotBFImages(mainP, data_DA, data_BF)

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
    thetas = mainP.getScanGrid(mainP.methods_set{m});
    radius = data_DA.Radius;
%     if mainP.interp_upsample > 0
%         radius = linspace(radius(1), radius(end), mainP.interp_upsample);
%     end
    warning('off')
    [scanConvertedImage, Xs, Zs] = getScanConvertedImage(m_BF, ...
        thetas, 1e3 * radius, 2024, 2024, 'spline');
    warning('on')
    img = db(abs(scanConvertedImage));
    xmin = Inf; xmax = -Inf; ymin = Inf; ymax = -Inf;
    for p=1:length(mainP.pts_range)
        ptraj = data_peaks{m}{p}.p_trajectory .* 1000;
        pxmin = min(ptraj(1,:)); pxmax = max(ptraj(1,:));
        pymin = min(ptraj(2,:)); pymax = max(ptraj(2,:));
        xmin = min(xmin, pxmin); xmax = max(xmax, pxmax);
        ymin = min(ymin, pymin); ymax = max(ymax, pymax);
    end

    clf(); subplot(1,2,1)
    imagesc(Xs, Zs, img)
    if mainP.normalize_bfim && (mainP.speckle_create || mainP.speckle_load)
        caxis([-25  25]);
    else
        caxis([max(img(:))-50 max(img(:))]);
    end
    xlim([xmin-1 xmax+1])
    ylim([ymin-1 ymax+1])
    colorbar
    colormap(gray)
    xlabel('azimuth [mm]');
    ylabel('range [mm]');
    title(['Method :', mainP.methods_set{m}])

    subplot(1,2,2); hold on
    img = db(abs(m_BF));
    pl_legend = {};
    angle_min = Inf; angle_max = -Inf; db_min = Inf; db_max = -Inf;
    for p=1:length(mainP.pts_range)
        pdata = data_peaks{m}{p};
        if p == 1
            pbtraj = data_peaks{m}{p}.beam_trajectory;
%             x_deg = rad2deg(pbtraj(1,:));
            x_mm = sin(pbtraj(1,:)) .* mainP.pts_range(p);
            plot(x_mm, pbtraj(3,:), colors_list{p}, 'LineWidth', 2, ...
                'LineStyle', linestyle_list{p}, 'Marker', markers_list{p});
            p_title = ['Steered response at ', ...
                num2str(mainP.pts_range(p), 2), ' mm range'];
            pl_legend{end+1} = p_title;
        end
        if isfield(pdata, 'peak') && isfield(pdata, 'peak_3db')
            p_y = pdata.peak(2) - 3;
%             p_x = rad2deg(pdata.peak_3db(1:2));
            p_x = sin(pdata.peak_3db(1:2)) .* mainP.pts_range(p);
            line('XData', p_x, 'YData', [p_y p_y], ...
                'LineWidth', 2, 'LineStyle', '-', 'Color', [0,0,0]+0.4)
%             p_l = strcat('P', int2str(p), '-3dB width= ', ...
%                 num2str(p_x(2) - p_x(1), 3), ' [deg]');
            p_l = strcat('P', int2str(p), '-3dB width= ', ...
                num2str(p_x(2) - p_x(1), 3), ' [mm]');
            pl_legend{end+1} = p_l;
            angle_min = min([angle_min, p_x(1)]);
            angle_max = max([angle_max, p_x(2)]);
            db_min = min([db_min, p_y]);
            db_max = max([db_max, p_y+3]);
        end
    end
    if angle_min < Inf && angle_max > -Inf
        xlim([angle_min-0.5, angle_max+0.5])
    end
    if db_min < Inf && db_max > -Inf
        ylim([db_min-3 db_max])
    end
%     xlabel('angle [deg]');
    xlabel('azimuth [mm]');
    ylabel('gain [dB]');
    legend(pl_legend, 'Location', 'south');
%     for t=1:length(mainP.P.Tx.Theta)
%         tx = rad2deg(mainP.P.Tx.Theta(t));
%         l = line('XData', [tx tx], 'YData', ylim, 'LineWidth', 2, ...
%             'LineStyle', linestyle_list{end-1}, ...
%             'Color', colors_list{end-1});
%     end
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
    thetas = mainP.getScanGrid(mainP.methods_set{m});
    radius = data_DA.Radius;
%     if mainP.interp_upsample > 0
%         radius = linspace(radius(1), radius(end), mainP.interp_upsample);
%     end
    warning('off')
%     [img, Xs, Zs] = getScanConvertedImage(m_BF, ...
%         thetas, 1e3 * radius, length(thetas), length(radius), 'spline'); % Test
    [img, Xs, Zs] = getScanConvertedImage(m_BF, ...
        thetas, 1e3 * radius, 2024, 2024, 'spline');
    warning('on')
    img = db(abs(img));
    
    shifts = mainP.shift.getShifts(mainP.P);
    shifts = interp1(1:length(shifts), shifts, ...
        linspace(1, length(shifts), length(thetas)));
    max_gain = max(img(:));
    az_lim = [NaN, -NaN];
    r_lim = [NaN, -NaN];
    for p=1:length(mainP.pts_range)
        for t=1:length(thetas)
            s_phantom = mainP.shift.shiftPositions(data_phantom, shifts(t));
            p_size = 0.1; %mm
            
            beam_az = atan(thetas(t)) * s_phantom.positions(p, 3)*1e3;
            pos_az = s_phantom.positions(p, 1)*1e3;
            pos_r = s_phantom.positions(p, 3)*1e3;
            s_az = find(abs(Xs - pos_az) < p_size, 10); % pos azimuth
%             s_az = find(abs(Xs - beam_az) < p_size, 10); % beam azimuth
            s_r = find(abs(Zs - pos_r) < p_size, 10);
            if ~isempty(s_r) && ~isempty(s_az) && abs(beam_az - pos_az) < 1
                % Create circle mask (instead of square)
%                 p_r = s_r * ones(1, length(s_az));
%                 p_az = s_az * ones(1, length(s_r));
%                 img_circle = zeros(size(img));
%                 p_circle = (p_r-mean(s_r)).^2 + (p_az'-mean(s_az)).^2 <= 10;
%                 img(s_r, s_az) = p_circle * (max_gain + 3) + ...
%                      abs(p_circle - 1).*img(s_r, s_az);
%                 img(s_r, s_az) = max_gain + 3; % square mask
                az_lim = [min(az_lim(1), pos_az - p_size), ...
                    max(az_lim(2), pos_az + p_size)];
                r_lim = [min(r_lim(1), pos_r - p_size), ...
                    max(r_lim(2), pos_r + p_size)];
            end
        end
    end
    pad = 0.5; % padding in mm
    az_lim = [max(Xs(1), az_lim(1) - pad), min(Xs(end), az_lim(2) + pad)];
    r_lim = [max(Zs(1), r_lim(1) - pad), min(Zs(end), r_lim(2) + pad)];
    minXs = find(Xs >= az_lim(1), 1);
    maxXs = find(Xs >= az_lim(2), 1);
    pXs = Xs(minXs:maxXs);
    minZs = find(Zs >= r_lim(1), 1);
    maxZs = find(Zs >= r_lim(end), 1);
    pZs = Zs(minZs:maxZs);
    p_img = img(minZs:maxZs, minXs:maxXs);
    p_img(~isfinite(p_img)) = -1000; % contourf doesn't handle non-finite values
%     [xGrid, zGrid] = meshgrid(linspace(pXs(1), pXs(end), 2024), ...
%         linspace(pZs(1), pZs(end), 2024));
    [xGrid, zGrid] = meshgrid(linspace(pXs(1), pXs(end), max([500, length(pXs)])), ...
        linspace(pZs(1), pZs(end), length(pZs)));
	p_img = interp2(pXs, pZs, p_img, xGrid, zGrid, 'spline', 0);
    subplot(2,ceil(length(mainP.methods_set)/2),m);
%     contourf(pXs, pZs, p_img, [-1000, max_gain-10, max_gain + 3], 'ShowText','on');
    contourf(xGrid, zGrid, p_img, [max_gain-100, max_gain-10, max_gain-3, max_gain + 3]) %, 'ShowText','on');
    colormap gray
    set(gca,'YDir','Reverse')
    xlabel('azimuth [mm]');
    ylabel('range [mm]');
    title(mainP.methods_set{m});
end
if mainP.save_plots
    saveas(gcf,  mainP.outputFileName('png'), 'png')
    saveas(gcf, mainP.outputFileName('fig'), 'fig')
    save(strcat(mainP.save_folder, mainP.files_prefix, ...
        'results_main_2_2.mat'), 'mainP', 'data_peaks', '-v7.3')
else
    pause
end
close
fprintf('Main_2_2 finished!\n')