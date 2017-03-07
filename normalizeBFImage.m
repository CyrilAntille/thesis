function [ bf_im ] = normalizeBFImage( m_peaks, bf_im, thetas, radius )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% background = anything outside beam_trajectory (+/- buffer)
traj_buffer = 0.5 * 1e-3; % m
bgd_mask = ones(size(bf_im));
bgd_mask(isnan(bf_im)) = 0;
for p=1:length(m_peaks)
    p_traj = m_peaks{p}.beam_trajectory; % [azimuth, radius, ampl]
    for t=1:size(p_traj, 1)
        min_angle = sin((p_traj(t,1) - traj_buffer)/p_traj(t,2));
        max_angle = sin((p_traj(t,1) + traj_buffer)/p_traj(t,2));
        angle_start = find(thetas >= min_angle, 1);
        angle_end = find(thetas >= max_angle, 1);
        radius_start = find(radius >= p_traj(t,2) - traj_buffer, 1);
        radius_end = find(radius >= p_traj(t,2) + traj_buffer, 1);
        bgd_mask(angle_start:angle_end, radius_start:radius_end) = 0;
    end
end
bgd_av = mean(bf_im(bgd_mask > 0));
bf_im = bf_im ./ bgd_av;

end

