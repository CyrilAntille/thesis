function [ bf_im ] = normalizeBFImage( mainP, bf_im, radius, phantom )
%normalizeBFImage Normalizes and/or upsamples BF images
%   Detailed explanation goes here
thetas = mainP.P.Tx.Theta;
if mainP.interp_upsample > 0
    warning('off')
    [thetaGrid, rangeGrid] = meshgrid(mainP.P.Tx.Theta, radius);
    thetas = linspace(mainP.P.Tx.Theta(1), mainP.P.Tx.Theta(end), mainP.interp_upsample);
%     radius = linspace(radius(1), radius(end), mainP.interp_upsample);
    % Comment out above line to upsample range dimension as well
    [upthetaGrid, uprangeGrid] = meshgrid(thetas, radius);
    bf_im = interp2(thetaGrid, rangeGrid, bf_im, ...
        upthetaGrid, uprangeGrid, mainP.interp_method, 0);
    warning('on')
end
if mainP.normalize_bfim
    scatterer_points = computePeaksInfo(mainP, phantom, radius, bf_im);
    % background = anything outside beam_trajectory (+/- buffer)
    traj_buffer = 0.5 * 1e-3; % m
    bgd_mask = ones(size(bf_im));
    bgd_mask(isnan(bf_im)) = 0;
    for p=1:length(scatterer_points)
        p_traj = scatterer_points{p}.beam_trajectory; % [azimuth, radius, ampl]
        for t=1:size(p_traj, 1)
            min_angle = asin((p_traj(t,1) - traj_buffer)/p_traj(t,2));
            max_angle = asin((p_traj(t,1) + traj_buffer)/p_traj(t,2));
            angle_start = find(thetas >= min_angle, 1);
            angle_end = find(thetas >= max_angle, 1);
            if isempty(angle_start)
                angle_start = 1;
            end
            if isempty(angle_end)
                angle_end = length(thetas);
            end
            radius_start = find(radius >= p_traj(t,2) - traj_buffer, 1);
            radius_end = find(radius >= p_traj(t,2) + traj_buffer, 1);
            bgd_mask(radius_start:radius_end, angle_start:angle_end) = 0;
        end
    end
    bgd_av = mean(bf_im(bgd_mask > 0));
    bf_im = bf_im ./ bgd_av;
%     bf_im = bf_im - bgd_av; % in dB -> substraction
end
% bf_im = db(abs(bf_im));
end

