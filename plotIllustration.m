function [] = plotIllustration( mainP, bfim, thetas, radius )

plot_velocities = false;

if mainP.save_plots
    figure('units','normalized','position',[.2 .3 .5 .3],'Visible','off')
%     set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 20*0.8 20*0.4])
else
    figure;
end

warning('off')
[scanConvertedImage, Xs, Zs] = getScanConvertedImage(bfim, ...
    thetas, 1e3 * radius, 2*2024, 2*2024, 'spline');
warning('on')
img = db(abs(scanConvertedImage));
% maxim = max(img(:)) + 100;
% color_range = [maxim-50, maxim];
maxim = 100;
color_range = [maxim-50, maxim];

pxl_size = 15;
% Transmitted beams
for t=1:length(mainP.P.Tx.Theta)
    x = Zs .* tan(mainP.P.Tx.Theta(t));
    for zidx=1:length(Zs)
        radius = sqrt(x(zidx)^2 + Zs(zidx)^2);
        [~, xidx] = min(abs(Xs - x(zidx)));
        if ~isempty(xidx) && xidx > 1 && xidx < length(Xs)
            if radius < mainP.medium_range(1) || ...
                    radius > mainP.medium_range(2) - 1
                continue % trims beams to image sector
            end
            xidx = xidx(1)-floor(pxl_size/2):xidx(1)+floor(pxl_size/2);
            xmin = max([xidx(1), 1]); xmax = min([xidx(end), length(Xs)]);
            xidx = xmin:xmax;
            img(zidx, xidx) = maxim;
        end
    end
end
% Image radius limits
lines_size = 0.05; % mm
for xidx=1:length(Xs)
    for zidx=1:length(Zs)
        radius = sqrt(Xs(xidx)^2 + Zs(zidx)^2);
        if abs(radius - mainP.medium_range(1)) < lines_size || ...
                abs(radius - mainP.medium_range(2) + 1) < lines_size
            img(zidx, xidx) = maxim;
        end
    end
end

if true %~plot_velocities
    % Points position
    p_size = 2.0; %mm
    disp_shift = Shift(ShiftType.RadialVar, 1/8, 17, 0, 1);
    shifts = disp_shift.getShifts(mainP.P)
    data_phantom = struct; data_phantom.positions = [0, 0, 40; 0, 0, 55]./1000;
    for s=1:length(shifts)
        s_phantom = mainP.shift.shiftPositions(data_phantom, shifts(s));
        for p=1:size(data_phantom.positions, 1)
            pos_az = s_phantom.positions(p, 1)*1e3;
            pos_r = s_phantom.positions(p, 3)*1e3;
            s_az = find(abs(Xs - pos_az) < p_size);
            s_r = find(abs(Zs - pos_r) < p_size);
            if ~isempty(s_r) && ~isempty(s_az)
                % img(s_r, s_az) = maxim; % square mask
                % Create circle mask (instead of square)
                p_r = s_r * ones(1, length(s_az));
                p_az = s_az * ones(1, length(s_r));
                p_circle = (p_r-mean(s_r)).^2 + (p_az'-mean(s_az)).^2 <= length(s_az);
                img(s_r, s_az) = p_circle * maxim + ...
                     abs(p_circle - 1).*img(s_r, s_az);
            end
        end
    end
end

% Plot image
clf(); imagesc(Xs, Zs, img)
xlabel('azimuth [mm]', 'FontSize', 14);
ylabel('range [mm]', 'FontSize', 14);
colormap(gray); caxis(color_range)

if plot_velocities
    hold on;
    % linestyle_list = {'-',':','-.','-.','-'}; % '   ', 
    linestyle_list = {'-','-','-','-','-','-','-','-','-'};
    % colors_list = {'b','r','g','m','c'};
    colors_list = {'w','w','w','w','w'};
    qlgd = {};
    
    % Random velocity vectors v_i
    norm_diag = norm([10, 10]);
    vs = [10, 0; [-10, 10]./norm_diag.*10; [-10, -10]./norm_diag.*10];
    vs_x = Xs(find(abs(Xs - 0) < 0.01, 1));
    vs_y = Zs(find(abs(Zs - 45) < 0.01, 1));
    for v=1:size(vs,1)
        quiver(vs_x, vs_y, vs(v,1), vs(v,2), 0.5, 'LineWidth', 5, ...
            'MaxHeadSize', 5, 'LineStyle', linestyle_list{v}, ...
            'Color', colors_list{v});
        vname = strcat('v_', int2str(v));
        text(vs_x+vs(v,1)/3, vs_y+vs(v,2)/3+1, vname, 'FontSize', 18, ...
            'FontWeight', 'bold', 'Color', 'w');
        qlgd{end+1} = horzcat(vname, ' = (', num2str(vs(v,1),2), ',', ...
            num2str(vs(v,2),2), ') m/s');
    end

    % Image acquisition lateral velocity v_tr
    vtr = [12, 0; 16.5, 0];
    vtr_range = [40, 55];
    for v=1:size(vtr,1)
        vtr_x = Xs(find(abs(Xs - 0) < 0.01, 1));
        vtr_y = Zs(find(abs(Zs - vtr_range(v)) < 0.01, 1));
        quiver(vtr_x, vtr_y, vtr(v,1), vtr(v,2), 0.5, 'LineWidth', 5, ...
            'MaxHeadSize', 5, 'LineStyle', linestyle_list{v}, ...
            'Color', colors_list{v});
        vname = strcat('v_{tr', int2str(v), '}');
        text(vtr_x+vtr(v,1)/3+1.5, vtr_y+vtr(v,2)/3+1, vname, 'FontSize', 18, ...
            'FontWeight', 'bold', 'Color', 'w');
        qlgd{end+1} = horzcat(vname, ' = (', num2str(vtr(v,1),3), ',', ...
            num2str(vtr(v,2),3), ') m/s');
    end

    % X-Y axes
    quiver(Xs(50), Zs(50), 5, 0, 'Color', 'w', 'LineWidth', 5)
    quiver(Xs(50), Zs(50), 0, 5, 'Color', 'w', 'LineWidth', 5)
    text(Xs(50)+5, Zs(50)+0.5, 'X', 'FontSize', 18, ...
        'FontWeight', 'bold', 'Color', 'w');
    text(Xs(50)+0.5, Zs(50)+5, 'Y', 'FontSize', 18,  ...
        'FontWeight', 'bold', 'Color', 'w');
    hold off; legend(qlgd, 'FontSize', 14);
end

if mainP.save_plots
    prefix = mainP.files_prefix;
    mainP.files_prefix = 'illustration_';
    saveas(gcf, mainP.outputFileName('png'), 'png')
    saveas(gcf, mainP.outputFileName('fig'), 'fig')
    mainP.files_prefix = prefix;
else
    pause
end
close

