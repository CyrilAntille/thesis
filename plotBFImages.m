function [] = plotBFImages( mainP, data_DA, data_BF )
%PLOTBFIMAGES Summary of this function goes here
%   Detailed explanation goes here
if mainP.save_plots
    figure('units','normalized','position',[.2 .3 .5 .3],'Visible','off')
else
    figure;
end

for m=1:length(mainP.methods_set)
    m_BF = data_BF{m};

    thetaRange = mainP.P.Tx.Theta;
    if strcmp(mainP.methods_set(m),'IAA-MBMB-Upsampled')
        thetaRange = linspace(m_P.Tx.Theta(1), m_P.Tx.Theta(end), ...
            mainP.upsample_number);
    end
    
    for s=1:mainP.shift.num_shifts
        if mainP.shift_per_beam
            s_DA = data_DA;
            s_BF = m_BF;
        else
            s_DA = data_DA{s};
            s_BF = m_BF{s};
        end
        warning('off')
        [scanConvertedImage, Xs, Zs] = getScanConvertedImage(s_BF, ...
            thetaRange, 1e3 * s_DA.Radius, 2024, 2024, 'spline');
        warning('on')
        img = db(abs(scanConvertedImage));
        clf(); imagesc(Xs, Zs, img)
        xlabel('azimuth [mm]');
    %             xlim([-15 15])
    %             ylim([35 45])
        if mainP.normalize_bfim
            caxis([-25  25]);
        end
        colorbar
        colormap(gray)
        ylabel('range [mm]');
        
        im_title = strcat('Method : ', mainP.methods_set{m});
        if mainP.shift_per_beam
            im_title = strcat(im_title, ', Shift: ', int2str(s));
            mainP.files_prefix = strcat('shift_', int2str(s), '_');
        end
        title(im_title)

        if mainP.save_plots
            prefix = mainP.files_prefix;
            mainP.files_prefix = strcat(mainP.files_prefix, ...
                mainP.methods_set{m}, '_', int2str(s), '_');
            output_file = mainP.outputFileName(true);
            saveas(gcf, output_file, 'png')
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

