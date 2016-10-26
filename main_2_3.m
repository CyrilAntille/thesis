%% 2.2: Motion between beams - 2 close points
% clear all
if ~exist('mainP', 'var')
    mainP = MainParameters();
    mainP.pts_range = [40 40]; % Add a range (in mm) for each point
    mainP.pts_azimuth = [0 4];
    mainP.num_beams = 101; % can be a single value or list of values
    mainP.shift = Shift(ShiftType.LinearSpeed, 0.5, -1, 0, 1);
%     mainP.shift = Shift(ShiftType.RadialVar, 7/8, -1, 0, 1);
    mainP.shift_per_beam = true;
    mainP.save_plots = false;
    
    mainP = mainP.createOutputDir();
end

%%
main_init

%% Dip measurements
dips = cell([length(mainP.num_beams), length(mainP.methods_set)]);
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
        % For this point on, Assumes 2 points at single range
        p_range = find(b_DA.Radius >= mainP.pts_range(1)*1e-3, 1);
        bf_img = bf_img(p_range,:);
%         [pdb, pdeg] = findpeaks(bf_img, thetaRange, 'SortStr','descend');
        [pdb, pidx] = findpeaks(bf_img, 'SortStr','descend');
        
        dip = struct; dip.depth = nan;
        if length(pdb) >= 2
            [pts_idx, idx] = sort([pidx(1) pidx(2)]);
            pts_gain = pdb(idx);
            dip.pts = [pts_idx; pts_gain];
            
            [min_gain, min_idx] = min(bf_img(pts_idx(1):pts_idx(2)));
            dip.dip_min = [pts_idx(1)+min_idx-1, min_gain];
            dip.depth = max(pts_gain) - min_gain;
            [min_peak, minpidx] = min(pts_gain);
            dip.scalloping = [pts_idx(minpidx) max(pts_gain) - min(pts_gain)];
        end
        dips{b,m} = dip;
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
        p_range = find(b_DA.Radius >= (mainP.pts_range(1))*1e-3, 1);
        bf_img = db(b_BF);
        bf_img = bf_img(p_range,:);
        thetaRange = rad2deg(thetaRange);
        azimuthRange = tand(thetaRange).*(mainP.pts_range(1));
        
        plot(azimuthRange, bf_img, 'LineWidth', 2);
        dip = dips{b,m};
        if ~isnan(dip.depth)
            dip_x = azimuthRange(dip.dip_min(1)); 
            line('XData', [dip_x dip_x], 'YData', [dip.dip_min(2) ...
                dip.dip_min(2)+dip.depth], 'LineWidth', 2, ...
                'LineStyle', linestyle_list{2}, 'Color', colors_list{2})
            p_l1 = strcat('Dip depth:  ', num2str(dip.depth, 3), ' dB');
        
            scallop_x = azimuthRange(dip.scalloping(1));
            line('XData', [scallop_x scallop_x], 'YData', [max(bf_img) ...
                - dip.scalloping(2) max(bf_img)], 'LineWidth', 2, ...
                'LineStyle', linestyle_list{1}, 'Color', colors_list{3})
            p_l2 = strcat('Scalloping loss:  ', ...
                num2str(dip.scalloping(2), 3), ' dB');
            
            pts_x = azimuthRange(dip.pts(1,:));
            line('XData', [pts_x(1), pts_x(2)], 'YData',...
                [max(dip.pts(2,:)) max(dip.pts(2,:))],  'LineWidth', 2, ...
                'LineStyle', linestyle_list{4}, 'Color', colors_list{4})
            p_l3 = strcat('Distance peaks:  ', ...
                num2str(pts_x(2) - pts_x(1), 3), ' degrees');
            
            legend({'Data', p_l1, p_l2, p_l3}, 'Location', 'best');
            ylim([dip.dip_min(2)-2 max(bf_img)+1])
%             xlim(pts_x + [-2, 2])
        end 
%         xlabel('angle [deg]');
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

