%% 2.2 - Motion between beams - Stats
clear all
mainP = MainParameters();
mainP.pts_range = [40, 50]; % Add a range (in mm) for each point
mainP.shift = Shift(ShiftType.RadialVar, 2/61, -1, 0); % Ref Shift.m
mainP.num_beams = 61; % can be a single value or list of values
mainP.shift_per_beam = true;

speeds = 0:0.5:3;

show_plots = true; % If false, saves plots to .png file
mainP = mainP.createOutputDir(~show_plots);

%%
pts_3dB_width = zeros([length(mainP.pts_range), length(mainP.num_beams),...
    length(mainP.methods_set), length(speeds)]);
for sp=1:length(speeds)
    mainP.shift = Shift(mainP.shift.type, speeds(sp), ...
        -1, mainP.shift.direction);
    main_init

    %3dB width computations
    max_peak_DAS = 0;
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
                pts_3dB_width(p,b,m,sp) = p_width;
            end
        end
    end
    clearvars -except mainP speeds pts_3dB_width
end

%% Plots
linestyle_list = {':','-','--','-.','-'};
markers_list = {'+','x','diamond','o','*'};

if show_plots
    figure;
else
    figure('units','normalized','position',[.2 .3 .5 .3],'Visible','off')
end

for p=1:size(pts_3dB_width, 1)
    for b=1:size(pts_3dB_width, 2)
        p1 = plot(speeds, squeeze(pts_3dB_width(p,b,:,:)), 'LineWidth', 2);
        for pidx=1:length(p1)
            p1(pidx).Marker = markers_list{pidx};
            p1(pidx).LineStyle = linestyle_list{pidx};
        end
        legend(mainP.methods_set, 'Location', 'best')
        title('Beampattern mainlobes 3dB width')
        xlabel('Speed [m/s]')
        ylabel('Mainlobe 3dB width [degrees]')
        if show_plots
            pause
        else
            im_name = strcat('speeds_p', int2str(mainP.pts_range(p)), '_',...
                int2str(mainP.num_beams(b)), '_', char(mainP.shift.type),...
                '_', int2str(mainP.shift.direction));
    %         saveas(gcf, strcat('../images/fig/', im_name, '.fig'), 'fig')
            saveas(gcf, strcat(mainP.save_folder, 'png/', im_name, '.png'), 'png')
        end
    end
end
close
fprintf('Stats_2_2 finished!')
