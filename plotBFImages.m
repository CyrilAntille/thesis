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
%         maxim = -50;
        color_range = [maxim-50, maxim];
        
        clf(); imagesc(Xs, Zs, img)
        xlabel('azimuth [mm]');
        if mainP.normalize_bfim && (mainP.speckle_create || mainP.speckle_load)
            color_range = [-25  25];
        end
        colorbar
        colormap(gray)
        ylabel('range [mm]');
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

