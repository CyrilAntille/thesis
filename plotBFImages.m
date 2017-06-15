function [] = plotBFImages( mainP, data_DA, data_BF, display_beams )
%PLOTBFIMAGES
if nargin < 4
    display_beams = false;
end

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
        
        if display_beams
%             maxim = maxim+100;
            for t=1:length(mainP.P.Tx.Theta)
                x = Zs .* tan(mainP.P.Tx.Theta(t));
                for zidx=1:length(Zs)
                    radius = sqrt(x(zidx)^2 + Zs(zidx)^2);
                    xidx = find(Xs >= x(zidx), 5);
                    if ~isempty(xidx) && radius >= mainP.medium_range(1) && ...
                            radius <= mainP.medium_range(2) - 1
                        if xidx(1) > ceil(length(xidx)/2)
                            xidx = xidx - ceil(length(xidx)/2);
                        end
                        img(zidx, xidx) = maxim;
                    end
                end
            end
            for xidx=1:length(Xs)
                for zidx=1:length(Zs)
                    radius = sqrt(Xs(xidx)^2 + Zs(zidx)^2);
                    if abs(radius - mainP.medium_range(1)) < 0.02 || ...
                            abs(radius - mainP.medium_range(2) + 1) < 0.02
                        img(zidx, xidx) = maxim;
                    end
                end
            end
            color_range = [maxim-50, maxim];
%             color_range = [-105  -55];
            
            % Display point position
            p_size = 2.0; %mm
            data_phantom = struct;
            data_phantom.positions = zeros([length(mainP.pts_range), 3]);
            shifts = mainP.shift.getShifts(mainP.P);
            if mainP.shift_per_beam
                shifts = interp1(1:length(shifts), shifts, ...
                    linspace(1, length(shifts), length(thetas)));
                for p=1:length(mainP.pts_range)
                    data_phantom.positions(p,:) = ...
                        [mainP.pts_azimuth(p), 0, mainP.pts_range(p)]./1000;
                    for t=1:length(shifts)
                        s_phantom = mainP.shift.shiftPositions(data_phantom, shifts(t));
                        pos_az = s_phantom.positions(p, 1)*1e3;
                        pos_r = s_phantom.positions(p, 3)*1e3;
                        s_az = find(abs(Xs - pos_az) < p_size, 1000);
                        s_r = find(abs(Zs - pos_r) < p_size, 1000);
                        if ~isempty(s_r) && ~isempty(s_az)
    %                         img(s_r, s_az) = maxim; % square mask
                            % Create circle mask (instead of square)
                            p_r = s_r * ones(1, length(s_az));
                            p_az = s_az * ones(1, length(s_r));
                            p_circle = (p_r-mean(s_r)).^2 + (p_az'-mean(s_az)).^2 <= length(s_az);
                            img(s_r, s_az) = p_circle * maxim + ...
                                 abs(p_circle - 1).*img(s_r, s_az);
                        end
                    end
                end
            else
                for p=1:length(mainP.pts_range)
                    data_phantom.positions(p,:) = ...
                        [mainP.pts_azimuth(p), 0, mainP.pts_range(p)]./1000;
                    pos_az = data_phantom.positions(p, 1)*1e3;
                    pos_r = data_phantom.positions(p, 3)*1e3;
                    s_az = find(abs(Xs - pos_az) < p_size, 1000);
                    s_r = find(abs(Zs - pos_r) < p_size, 1000);
                    if ~isempty(s_r) && ~isempty(s_az)
%                         img(s_r, s_az) = maxim; % square mask
                        % Create circle mask (instead of square)
                        p_r = s_r * ones(1, length(s_az));
                        p_az = s_az * ones(1, length(s_r));
                        p_circle = (p_r-mean(s_r)).^2 + (p_az'-mean(s_az)).^2 <= length(s_az);
                        img(s_r, s_az) = p_circle * maxim + ...
                             abs(p_circle - 1).*img(s_r, s_az);
                    end
                end
            end
            % end display position
        end
        
        clf(); imagesc(Xs, Zs, img)
        xlabel('azimuth [mm]');
        if mainP.normalize_bfim && (mainP.speckle_create || mainP.speckle_load)
            color_range = [-25  25];
        end
        colorbar
        colormap(gray)
        ylabel('range [mm]');
        caxis(color_range)
        
        if display_beams
            linestyle_list = {'-',':','-.','-.','-'};
            colors_list = {'b','r','g','m','c'};
            norm_diag = norm([0.5, 0.5]);
            vs = [0.5, 0; [-0.5, 0.5]./norm_diag./2; [-0.5, -0.5]./norm_diag./2];
            vs_x = Xs(find(abs(Xs - 0) < 0.01, 1)); 
            vs_y = Zs(find(abs(Zs - 45) < 0.01, 1));
            hold on;
            qlgd = cell([1 size(vs,1)+2]);
            for v=1:size(vs,1)
                quiver(vs_x, vs_y, vs(v,1), vs(v,2), 4, 'LineWidth', 5, ...
                    'MaxHeadSize', 5, 'LineStyle', linestyle_list{v}, ...
                    'Color', colors_list{v});
                qlgd{v} = horzcat('   v_s = (', num2str(vs(v,1),2), ',', ...
                    num2str(vs(v,2),2), ') m/s');
            end

            vtr_40 = [12, 0];
            vtr_40_y = Zs(find(abs(Zs - 40) < 0.01, 1));
            quiver(vs_x, vtr_40_y, vtr_40(1), vtr_40(2), 0.5, 'LineWidth', 5, ...
                'MaxHeadSize', 5, 'LineStyle', linestyle_list{v+1}, ...
                'Color', colors_list{v+1});
            qlgd{v+1} = horzcat('   v_{tr} = (', num2str(vtr_40(1),3), ',', ...
                num2str(vtr_40(2),2), ') m/s');

            vtr_55 = [16.5, 0];
            vtr_55_y = Zs(find(abs(Zs - 55) < 0.01, 1));
            quiver(vs_x, vtr_55_y, vtr_55(1), vtr_55(2), 0.5, 'LineWidth', 5, ...
                'MaxHeadSize', 5, 'LineStyle', linestyle_list{v+2}, ...
                'Color', colors_list{v+2});
            qlgd{v+2} = horzcat('   v_{tr} = (', num2str(vtr_55(1),3), ',', ...
                num2str(vtr_55(2),2), ') m/s');

            quiver(Xs(50), Zs(50), 5, 0, 'Color', 'w', 'LineWidth', 5)
            quiver(Xs(50), Zs(50), 0, 5, 'Color', 'w', 'LineWidth', 5)
    %         qlgd{end+1} = 'X'; qlgd{end+1} = 'Y';
            text(Xs(50)+5, Zs(50)+0.5, 'X', 'FontSize', 18, ...
                'FontWeight', 'bold', 'Color', 'w');
            text(Xs(50)+0.5, Zs(50)+5, 'Y', 'FontSize', 18,  ...
                'FontWeight', 'bold', 'Color', 'w');
            hold off; legend(qlgd, 'FontSize', 12);
        end
        
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

