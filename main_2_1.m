%% 2.1: Motion between frames
clear all
mainP = MainParameters();
mainP.pts_range = [40, 50]; % Add a range (in mm) for each point
mainP.shift = Shift(ShiftType.LinearCst, 2*1e-3, 3, 0); % Ref Shift.m
mainP.num_beams = 101; % can be a single value or list of values
mainP.shift_per_beam = false;

mainP.save_plots = false;
grayscale_plots = false;
mainP = mainP.createOutputDir();

%%
main_init

%% Estimate scalloping loss
fprintf('Scalloping loss estimation and plotting\n')
theta_max = asin(mainP.P.Tx.SinThMax);
% Maximum scalloping loss (for all beams)
loss_beams = zeros([length(mainP.pts_range) ...
    length(mainP.methods_set) length(mainP.num_beams)]);
% Scalloping loss per shift (for a given number of beams)
chosen_num_beams = mainP.num_beams(1);
loss_shift = zeros([length(mainP.pts_range) ...
    length(mainP.methods_set) length(mainP.shift.num_shifts)]);

for m=1:length(mainP.methods_set)
    m_BF = data_BF{m};
    fprintf('----------------------\n')
    fprintf('Beamforming method: %s.\n', mainP.methods_set{m})
    for b=1:length(mainP.num_beams)
        b_BF = m_BF{b};
        b_DA = data_DA{b};
        pts_ampl = ones(mainP.shift.num_shifts,length(mainP.pts_range))*Inf;
        for s=1:mainP.shift.num_shifts
            s_BF = b_BF{s};
            s_DA = b_DA{s};
            speckle_pts = data_phantoms{b}{s};
            [Theta,Phi,R] = cart2sph(speckle_pts.positions(:,3),...
                speckle_pts.positions(:,1),speckle_pts.positions(:,2));
            
            radiusRange = s_DA.Radius;
            thetaRange = linspace(-theta_max, theta_max, mainP.num_beams(b));
            if strcmp(mainP.methods_set(m),'IAA-MBMB-Upsampled')
                thetaRange = linspace(-theta_max, theta_max, mainP.upsample_number);
            end
            s_pts = pts_ampl(s,:);
            for p=1:length(Theta)
                pt = Theta(p);
                pr = R(p);
                s_pts(p) = Inf;
                if pt < thetaRange(1) || pt > thetaRange(end) || ...
                        pr < radiusRange(1) || pr > radiusRange(end)
                    continue;
                end
                [~, tidx] = min(abs(pt-thetaRange));
                [~, ridx] = min(abs(pr-radiusRange));
                s_pts(p) = db(abs(s_BF(ridx,tidx)));
            end
            pts_ampl(s,:) = s_pts;
            if mainP.num_beams(b) == chosen_num_beams
                loss_shift(:,m,s) = s_pts;
            end
            
        end
        for p=1:size(pts_ampl, 2)
            if max(pts_ampl(:,p)) == Inf
                continue % -> point not in beamformed area
            end
            mloss = max(pts_ampl(:,p)) - min(pts_ampl(:,p));
            loss_beams(p,m,b) = mloss;
        end
    end
end

%% Plots
linestyle_list = {':','-','--','-.','-'};
markers_list = {'+','x','diamond','o','*'};
line_clr = 'green';
if grayscale_plots
    line_clr = [0,0,0]+0.6;
end

if mainP.save_plots
    figure('units','normalized','position',[.2 .3 .5 .3],'Visible','off')
else
    figure;
end
% Loss vs shift
for p=1:length(mainP.pts_range)
    shifts = (0:mainP.shift.num_shifts-1) * mainP.shift.val;
    pl = plot(shifts, squeeze(loss_shift(p,:,:)), 'LineWidth', 2);
    grayscale_alpha = linspace(0,0.4,length(pl));
    for pidx=1:length(pl)
        pl(pidx).Marker = markers_list{pidx};
        pl(pidx).LineStyle = linestyle_list{pidx};
        if grayscale_plots
            pl(pidx).Color = [0,0,0]+grayscale_alpha(pidx);
        end
    end
    b = 0;
    while true
        if b - 1 > max(shifts)
            break
        end
        l = line('XData', [b b], 'YData', ylim, ...
            'LineWidth', 2, 'LineStyle', '-.', 'Color', line_clr);
        b = b + 1;
    end
    legend([mainP.methods_set, 'Transmitted beams'], 'Location', 'best');
    ylabel('Scalloping loss [dB]');
    xlabel('Shift [ratio beams separation]');
    t = strcat('Scatterer point at ', num2str(mainP.pts_range(p),0), 'mm range, ');
%     title(t)

    if mainP.save_plots
        im_name = strcat('loss_shift_p', int2str(mainP.pts_range(p)), ...
            '_', int2str(chosen_num_beams), '_', char(mainP.shift.type),...
            '_', int2str(mainP.shift.direction));
        saveas(gcf, strcat(mainP.save_folder, 'png/', im_name, '.png'), 'png')
    else
        pause
    end
end

% Max loss vs beams
for p=1:length(mainP.pts_range)
    pl = plot(mainP.num_beams, squeeze(loss_beams(p,:,:))', 'LineWidth', 2);
    grayscale_alpha = linspace(0,0.4,length(pl));
    for pidx=1:length(pl)
        pl(pidx).Marker = markers_list{pidx};
        pl(pidx).LineStyle = linestyle_list{pidx};
        if grayscale_plots
            pl(pidx).Color = [0,0,0]+grayscale_alpha(pidx);
        end
    end

    line('XData', [mainP.num_beams(1) mainP.num_beams(end)], 'YData', ...
        [1 1], 'LineWidth', 2, 'LineStyle', '-.', 'Color', line_clr);
    legend([mainP.methods_set '1dB threshold'], 'Location', 'best');
    ylabel('Max scalloping loss [dB]');
    xlabel('Number of transmitted beams');
    if length(mainP.num_beams) > 1
        xlim([mainP.num_beams(1) mainP.num_beams(end)]) 
    end
%     ylim([0 5])
    if mainP.save_plots
        im_name = strcat('loss_beams_p', int2str(mainP.pts_range(p)), ...
            '_', char(mainP.shift.type), '_', int2str(mainP.shift.direction));
        saveas(gcf, strcat(mainP.save_folder, 'png/', im_name, '.png'), 'png')
    else
        pause
    end
end

%% Beamformed images plots
if ~mainP.save_plots
    figure;
    for m=1:length(mainP.methods_set)
        m_BF = data_BF{m};
        for b=1:length(mainP.num_beams)
            Pb = mainP.copyP(mainP.num_beams(b));
            b_BF = m_BF{b};
            b_DA = data_DA{b};
            for s=1:mainP.shift.num_shifts
                s_BF = b_BF{s};
                s_DA = b_DA{s};

                thetaRange = Pb.Tx.Theta;
                if strcmp(mainP.methods_set(m),'IAA-MBMB-Upsampled')
                    thetaRange = linspace(Pb.Tx.Theta(1), Pb.Tx.Theta(end), ...
                        mainP.upsample_number);
                end
                warning('off')
                [scanConvertedImage, Xs, Zs] = getScanConvertedImage(s_BF, ...
                    thetaRange, 1e3 * s_DA.Radius, 2024, 2024, 'spline');
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
                    int2str(mainP.num_beams(b)), ', Shift No: ', int2str(s-1)])
                pause
            end
        end
    end
    close
end
fprintf('Main_2_1 finished!')