function [] = plotBFImages( mainP, data_DA, data_BF )
%PLOTBFIMAGES
if mainP.save_plots
    figure('units','normalized','position',[.2 .3 .5 .3],'Visible','off')
else
    figure;
end
for m=1:length(mainP.methods_set)
    thetas = mainP.getScanGrid(mainP.methods_set{m});
    m_BF = data_BF{m};
    for s=1:mainP.shift.num_shifts
        if mainP.shift_per_beam
            s_BF = m_BF;
            radius = data_DA.Radius;
        else
            s_BF = m_BF{s};
            radius = data_DA{s}.Radius;
        end
        if mainP.interp_upsample > 0
            radius = linspace(radius(1), radius(end), mainP.interp_upsample);
        end
        
        warning('off')
        [scanConvertedImage, Xs, Zs] = getScanConvertedImage(s_BF, ...
            thetas, 1e3 * radius, 2024, 2024, 'spline');
        warning('on')
        img = db(abs(scanConvertedImage));
        maxim = max(img(:));
        maxim = -50;
        color_range = [maxim-50, maxim];
        
        
        % Points position
        p_size = 2.0; %mm
%         disp_shift = Shift(ShiftType.RadialVar, 1/8, 17, 0, 1);
        shifts = mainP.shift.getShifts(mainP.P);
        data_phantom = struct; data_phantom.positions = [0, 0, 40; 0.75, 0, 40]./1000;
        for sh=1:length(shifts)
            s_phantom = mainP.shift.shiftPositions(data_phantom, shifts(sh));
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
        
        
        
        clf(); imagesc(Xs, Zs, img)
        xlabel('azimuth [mm]', 'FontSize', 14);
        if mainP.normalize_bfim && (mainP.speckle_create || mainP.speckle_load)
            color_range = [-25  25];
        end
        colorbar
        colormap(gray)
        ylabel('range [mm]', 'FontSize', 14);
        caxis(color_range)
        
        im_title = strcat('Method : ', mainP.methods_set{m});
        if mainP.shift_per_beam
            im_title = strcat(im_title, ', Shift: ', int2str(s));
            mainP.files_prefix = strcat('shift_', int2str(s), '_');
        end
%         title(im_title)

        if mainP.save_plots
            prefix = mainP.files_prefix;
            mainP.files_prefix = strcat(mainP.files_prefix, ...
                mainP.methods_set{m}, '_', int2str(s), '_');
            saveas(gcf, mainP.outputFileName('png'), 'png')
            saveas(gcf, mainP.outputFileName('fig'), 'fig')
            mainP.files_prefix = prefix;
        else
            pause
        end
        if mainP.shift_per_beam
            break
        end
    end
end
close

