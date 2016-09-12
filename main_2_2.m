%% 2.2: Motion between beams
clear all
mainP = MainParameters();
mainP.pts_range = [40, 50]; % Add a range (in mm) for each point
mainP.num_beams = 101; % can be a single value or list of values
mainP.shift = Shift(ShiftType.LinearSpeed, 1, -1, -90);
% mainP.shift = Shift(ShiftType.RadialVar, 1/2, -1, 0);
mainP.shift_per_beam = true;

show_plots = false; % If false, saves plots to .png file
grayscale_plots = false;
mainP = mainP.createOutputDir(~show_plots);

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
    points_speed = (shifts_p(end) - shifts_p(1)) / time_frame % m/s
end

%%
main_init

%% 3dB width computations

max_peak_DAS = 0;
pts_3dB = cell([length(mainP.pts_range), ...
    length(mainP.methods_set), length(mainP.num_beams)]);
inters_3dB = cell([length(mainP.pts_range), ...
    length(mainP.methods_set), length(mainP.num_beams)]);
for m=1:length(mainP.methods_set)
    m_BF = data_BF{m};
    for b=1:length(mainP.num_beams)
        Pb = mainP.copyP(mainP.num_beams(b));
        b_BF = m_BF{b};
        b_DA = data_DA{b};
        if m == 1
            max_peak_DAS = max(db(abs(b_BF(:))));
        end

        thetaRange = Pb.Tx.Theta;
        if strcmp(mainP.methods_set(m),'IAA-MBMB-Upsampled')
            thetaRange = linspace(Pb.Tx.Theta(1), Pb.Tx.Theta(end), 500);
        end
        
        thetaRange = rad2deg(thetaRange);
        bf_img = db(b_BF)  - max_peak_DAS;
        for p=1:length(mainP.pts_range)
            p_range = mainP.pts_range(p); % range not constant if LateralShift
            minr = find(b_DA.Radius >= (p_range - 2)*1e-3, 1);
            maxr = find(b_DA.Radius >= (p_range + 2)*1e-3, 1);
            p_bp = max(bf_img(minr:maxr, :), [], 1);
            [pdb, pdeg] = findpeaks(p_bp, thetaRange, ...
                'SortStr','descend');
            p_3dbline = ones(size(thetaRange)) .* (pdb(1) - 3);
            [p_x, p_y, ~, ~] = intersections(thetaRange, p_bp, ...
                thetaRange, p_3dbline, 1);
            if length(p_x) < 2
                p_width = NaN;
            else
                p_width = p_x(2) - p_x(1);
            end
            pts_3dB{p,m,b} = [pdb(1), p_width];
            inters_3dB{p,m,b} = p_x;
        end
    end
end
        

%% Plots
linestyle_list = {'-.','--','-',':'};
markers_list = {'+','x','diamond','o'};
colors_list = {'b','r','g','k','m','y'};

if show_plots
    figure;
else
    figure('units','normalized','position',[.2 .3 .5 .3],'Visible','off')
end
max_peak_DAS = 0; max_peak_DAS_img = 0;
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
        
        if m == 1
            max_peak_DAS = max(db(abs(b_BF(:))));
            max_peak_DAS_img = max(img(:));
        end
        
        cla(); subplot(1,2,1)
        imagesc(Xs, Zs, img - max_peak_DAS_img)
        caxis([-50  0]);
%         imagesc(Xs, Zs, img)
%         caxis([-120  -70]);
        xlim([-20 20])
        ylim([min(mainP.pts_range)-5 max(mainP.pts_range)+5])
        colorbar
        colormap(gray)
        xlabel('azimuth [mm]');
        ylabel('range [mm]');
        title(['Method :', mainP.methods_set{m}, ', No Beams: ', ...
            int2str(mainP.num_beams(b))])

        
        subplot(1,2,2); hold on
        bf_img = db(b_BF) - max_peak_DAS;
        thetaRange = rad2deg(thetaRange);
        pl_legend = {};
        y_min = Inf; x_min = Inf; x_max = -Inf;
        for p=1:length(mainP.pts_range)
            p_range = mainP.pts_range(p); % range not constant if LateralShift
            minr = find(b_DA.Radius >= (p_range - 2)*1e-3, 1);
            maxr = find(b_DA.Radius >= (p_range + 2)*1e-3, 1);
            p_bp = max(bf_img(minr:maxr, :), [], 1);
            plot(thetaRange, p_bp, colors_list{p}, 'LineWidth', 2, ...
                'LineStyle', linestyle_list{p}, 'Marker', markers_list{p});
            
            p_y = pts_3dB{p,m,b}(1) - 3;
            p_title = strcat('P', int2str(p), '-', ...
                int2str(data_phantoms.positions(p,3)*1000), 'mm range');
            line('XData', inters_3dB{p,m,b}, 'YData', [p_y p_y], ...
                'LineWidth', 2, 'LineStyle', '-', 'Color', [0,0,0]+0.4)
            p_l = strcat('P', int2str(p), '-3dB width= ', ...
                num2str(pts_3dB{p,m,b}(2), 3), ' degrees');
            
            pl_legend{end+1} = p_title;
            pl_legend{end+1} = p_l;
            x_min = min([x_min, inters_3dB{p,m,b}(1)]);
            x_max = max([x_max, inters_3dB{p,m,b}(2)]);
            y_min = min([y_min, p_y]);
        end
        xlim([x_min-1, x_max+1])
        ylim([y_min-3 0])
        xlabel('angle [deg]');
        ylabel('gain [dB]');
        legend(pl_legend, 'Location', 'best');
        hold off;
        
        if show_plots
            pause
        else
            im_name = strcat(int2str(mainP.num_beams(b)), '_', ...
                char(mainP.shift.type), '_', int2str(mainP.shift.direction),...
                '_', mainP.methods_set{m});
    %         saveas(gcf, strcat('../images/fig/', im_name, '.fig'), 'fig')
            saveas(gcf, strcat(mainP.save_folder, 'png/', im_name, '.png'), 'png')
        end
    end
end
close
fprintf('\nMain_2_2 finished!\n')
%% Beamformed images plots
if show_plots
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
            caxis([-120  -70]);
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