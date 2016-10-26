%% 2.2: Motion between beams
% clear all
if ~exist('mainP', 'var')
    mainP = MainParameters();
    mainP.pts_range = [40 50]; % Add a range (in mm) for each point
    mainP.pts_azimuth = [0 1];
    mainP.num_beams = 101; % can be a single value or list of values
    mainP.shift = Shift(ShiftType.LinearSpeed, 1, -1, 90, 1);
%     mainP.shift = Shift(ShiftType.RadialVar, 7/8, -1, 0, 1);
    mainP.shift_per_beam = true;
    mainP.save_plots = false;
    
    mainP = mainP.createOutputDir();
end

%% Test - scatterer points speed estimation in m/s
time_frame = 0.2 * 1e-3 * mainP.num_beams(1); % s
Ptest = mainP.copyP(mainP.num_beams(1));
b_shift = Shift(mainP.shift.type, mainP.shift.val, ...
                Ptest.Tx.NTheta, mainP.shift.direction);
shifts = b_shift.getShifts(Ptest); % in degrees or m
for p=1:length(mainP.pts_range)
    if mainP.shift.type == ShiftType.RadialVar || ...
            mainP.shift.type == ShiftType.RadialCst
        shifts_p = sin(shifts) * mainP.pts_range(p) * 1e-3; % m (approx)
    else
        shifts_p = shifts; % m
    end
    points_speed = (shifts_p(end) - shifts_p(1)) / time_frame; % m/s
end

%%
main_init

%% 3dB width computations
pts_3dB = cell([length(mainP.pts_range), ...
    length(mainP.num_beams), length(mainP.methods_set)]);
inters_3dB = cell([length(mainP.pts_range), ...
    length(mainP.num_beams), length(mainP.methods_set)]);
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
        for p=1:length(mainP.pts_range)
            p_range = mainP.pts_range(p); % range not constant if LateralShift
            minr = find(b_DA.Radius >= (p_range - 2)*1e-3, 1);
            maxr = find(b_DA.Radius >= (p_range + 2)*1e-3, 1);
            p_bp = max(bf_img(minr:maxr, :), [], 1);
            [pdb, pdeg] = findpeaks(p_bp, thetaRange, 'SortStr','descend');
            p_3dbline = ones(size(thetaRange)) .* (pdb(1) - 3);
            [p_x, p_y, ~, ~] = intersections(thetaRange, p_bp, ...
                thetaRange, p_3dbline, 1);
            if length(p_x) < 2
                p_width = NaN;
            elseif length(p_x) > 2
                p_x2 = max(find(p_x >= pdeg(1), 1), 2);
                p_width = p_x(p_x2) - p_x(p_x2-1);
                p_x = p_x(p_x2-1:p_x2);
            else
                p_width = p_x(2) - p_x(1);
            end
            pts_3dB{p,b,m} = [pdb(1), p_width];
            inters_3dB{p,b,m} = p_x;
        end
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
        bf_img = db(abs(b_BF));
        thetaRange = rad2deg(thetaRange);
        pl_legend = {};
        y_min = Inf; y_max = -Inf; x_min = Inf; x_max = -Inf;
        for p=1:length(mainP.pts_range)
            p_range = mainP.pts_range(p); % range not constant if LateralShift
            minr = find(b_DA.Radius >= (p_range - 2)*1e-3, 1);
            maxr = find(b_DA.Radius >= (p_range + 2)*1e-3, 1);
            p_bp = max(bf_img(minr:maxr, :), [], 1);
            plot(thetaRange, p_bp, colors_list{p}, 'LineWidth', 2, ...
                'LineStyle', linestyle_list{p}, 'Marker', markers_list{p});
            
            p_y = pts_3dB{p,b,m}(1) - 3;
            p_title = strcat('P', int2str(p), '-', ...
                int2str(data_phantoms.positions(p,3)*1000), 'mm range');
            pl_legend{end+1} = p_title;
            if length(inters_3dB{p,b,m}) == 2
                line('XData', inters_3dB{p,b,m}, 'YData', [p_y p_y], ...
                    'LineWidth', 2, 'LineStyle', '-', 'Color', [0,0,0]+0.4)
                p_l = strcat('P', int2str(p), '-3dB width= ', ...
                    num2str(pts_3dB{p,b,m}(2), 3), ' degrees');

                pl_legend{end+1} = p_l;
                x_min = min([x_min, inters_3dB{p,b,m}(1)]);
                x_max = max([x_max, inters_3dB{p,b,m}(2)]);
                y_min = min([y_min, p_y]);
                y_max = max([y_max, pts_3dB{p,b,m}(1)]);
            end
        end
        xlim([x_min-1, x_max+1])
        ylim([y_min-3 y_max])
        xlabel('angle [deg]');
        ylabel('gain [dB]');
        legend(pl_legend, 'Location', 'best');
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

for b=1:length(mainP.num_beams)
    for m=1:length(mainP.methods_set)
        Pb = mainP.copyP(mainP.num_beams(b));
        b_BF = data_BF{m}{b};
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
        
        
        if mainP.shift_per_beam
            b_shift = Shift(mainP.shift.type, mainP.shift.val, ...
                Pb.Tx.NTheta, mainP.shift.direction);
            beam_shift = Pb.Tx.SinTheta(2) - Pb.Tx.SinTheta(1);
        else
            b_shift = mainP.shift;
        end
        shifts = b_shift.getShifts(Pb);
        max_gain = max(img(:));
        for p=1:length(mainP.pts_range)
            for s=1:b_shift.num_shifts
%             for s=floor(b_shift.num_shifts*1/3):ceil(b_shift.num_shifts*2/3)
                s_phantom = b_shift.shiftPositions(data_phantoms, shifts(s));
                p_size = 0.05; %mm
%                 s_az = find(abs(Xs - s_phantom.positions(p, 1)*1e3) < p_size, 10);
                beam_az = atan(Pb.Tx.Theta(s)) * s_phantom.positions(p, 3)*1e3;
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
            contourf(pXs, pZs, p_img, [-100, pts_3dB{p,b,m}(1)-10, ...
                pts_3dB{p,b,m}(1)-3, max_gain + 3], 'ShowText','on');
            set(gca,'YDir','Reverse')
            xlabel('azimuth [mm]');
            ylabel('range [mm]');
            title(mainP.methods_set{m});
        end
    end
    if mainP.save_plots
        im_name = strcat(int2str(mainP.num_beams(b)), '_', ...
            char(mainP.shift.type), '_', int2str(mainP.shift.direction),...
            '_', num2str(mainP.shift.val,2), '_p');
        for p=1:length(mainP.pts_range)
            p_name = strcat(im_name, int2str(mainP.pts_range(p)), 'mm');
            saveas(figs(p), strcat(mainP.save_folder, 'png/', p_name, '.png'), 'png')
        end
    else
        pause
    end
end
close all
set(0, 'DefaultFigureVisible', 'on')

%% Beamformed images plots
if ~mainP.save_plots
    figure;
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
            imagesc(Xs, Zs, img)
            xlabel('azimuth [mm]');
    %             xlim([-15 15])
    %             ylim([35 45])
            caxis([-25  25]);
            colorbar
            colormap(gray)
            ylabel('range [mm]');
            title(['Method :', mainP.methods_set{m}, ', No Beams: ', ...
                int2str(mainP.num_beams(b))])
            pause
        end
    end
    close
end
fprintf('Main_2_2 finished!\n')