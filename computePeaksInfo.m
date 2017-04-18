function [ scatterer_points ] = computePeaksInfo(mainP, phantom, da_radius, bf_im)
%COMPUTEPEAKSINFO Summary of this function goes here
%   Detailed explanation goes here

bf_im = db(abs(bf_im));
scatterer_points = cell([1 length(mainP.pts_range)]);
pts_center = zeros([1 length(mainP.pts_range)]);
shifts = mainP.shift.getShifts(mainP.P);

% Fits radius and thetas to image size
radius = linspace(da_radius(1), da_radius(end), size(bf_im, 1));
shifts = linspace(shifts(1), shifts(end), size(bf_im, 2));
thetas = linspace(mainP.P.Tx.Theta(1), mainP.P.Tx.Theta(end), size(bf_im, 2));
for p=1:length(mainP.pts_range)
    scat_p = struct;
    scat_p.p_trajectory = zeros([3 length(shifts)]); % [x, z, ampl]
    scat_p.beam_trajectory = zeros([3 length(shifts)]); % [angle, radius, ampl]
    % Scatterer point and beam trajectories
    if mainP.shift_per_beam
        for s=1:length(shifts)
            s_phantom = mainP.shift.shiftPositions(phantom, shifts(s));
            scat_p.p_trajectory(:,s) = [s_phantom.positions(p, 1), ...
                s_phantom.positions(p, 3), s_phantom.amplitudes(p)];

            s_az = tan(thetas(s)) * s_phantom.positions(p, 3);
            s_radius = sqrt(s_phantom.positions(p, 3)^2 + s_az^2);
            s_radius_idx = find(radius >= s_radius - 5 * 1e-4, 1);
            s_radius_end = find(radius >= s_radius + 5 * 1e-4, 1);
            if ~isempty(s_radius_idx) && isempty(s_radius_end)
                s_radius_end = length(radius);
            elseif isempty(s_radius_idx) && ~isempty(s_radius_end)
                s_radius_idx = 1;
            end
            if ~isempty(s_radius_idx)
                bf_ampl = max(bf_im(s_radius_idx:s_radius_end, s));
%                 scat_p.beam_trajectory(:,s) = [s_az; s_radius; bf_ampl];
                scat_p.beam_trajectory(:,s) = [thetas(s); s_radius; bf_ampl];
            end
        end
    else
        point_pos = [phantom.positions(p, 1); ...
            phantom.positions(p, 3); phantom.amplitudes(p)];
        scat_p.p_trajectory = repmat(point_pos, 1, length(thetas));
        s_radius = sqrt(phantom.positions(p, 1)^2 + ...
            phantom.positions(p, 3)^2);
        s_radius_idx = find(radius >= s_radius - 5 * 1e-4, 1);
        s_radius_end = find(radius >= s_radius + 5 * 1e-4, 1);
        if isempty(s_radius_end)
            s_radius_end = length(radius);
        end
        if ~isempty(s_radius_idx)
            scat_p.beam_trajectory = ones(size(scat_p.beam_trajectory)) .* s_radius;
%             scat_p.beam_trajectory(1,:) = tan(thetas) .* phantom.positions(p, 3);
            scat_p.beam_trajectory(1,:) = tan(thetas);
            scat_p.beam_trajectory(3,:) =  ...
                max(bf_im(s_radius_idx:s_radius_end, :), [], 1);
        end
    end
    % Expected point angle when beam hits it
    p_angle = asin(scat_p.p_trajectory(1,:) ./ scat_p.p_trajectory(2,:));
    diff_angle = abs(p_angle - scat_p.beam_trajectory(1,:));
    [~, angle_idx] = min(diff_angle);
    pts_center(p) = p_angle(angle_idx);
    scatterer_points{p} = scat_p;
end

for p=1:length(mainP.pts_range)
    scat_p = scatterer_points{p};
    % Finds nearest peak
    [pdb, pangle] = findpeaks(scat_p.beam_trajectory(3,:), scat_p.beam_trajectory(1,:));
    if isempty(pangle)
        break
    end
    [~, pangle_idx] = min(abs(pangle - pts_center(p)));
    [~, pt_idx] = min(abs(pts_center - pangle(pangle_idx)));
    if ~mainP.shift_per_beam || (pts_center(p) - pangle(pangle_idx) < 1e-3 && (pt_idx == p))
        if pangle(pangle_idx) - pts_center(p) < 0 && pangle_idx < length(pangle)
            % -> compare with next index
            [~, pt_idx] = min(abs(pts_center - pangle(pangle_idx+1)));
            if pdb(pangle_idx+1) > pdb(pangle_idx) && pt_idx == p && ...
                    pts_center(p) - pangle(pangle_idx + 1) < 25*1e-3
                pangle_idx = pangle_idx + 1; % 25*1e-3 rad = atan(1/40 mm)
            end
        elseif pangle(pangle_idx) - pts_center(p) >= 0 && pangle_idx > 1
            % -> compare with previous index
            [~, pt_idx] = min(abs(pts_center - pangle(pangle_idx-1)));
            if pdb(pangle_idx-1) > pdb(pangle_idx) && pt_idx == p && ...
                   pts_center(p) - pangle(pangle_idx - 1) < 25*1e-3
                pangle_idx = pangle_idx - 1;
            end
        end
        scat_p.peak = [pangle(pangle_idx); pdb(pangle_idx)]; % [angle, ampl]
        % Peak 3dB width
        p_3dbline = ones([1 size(scat_p.beam_trajectory,2)]).*(scat_p.peak(2) - 3);
        [p_cross, ~, ~, ~] = intersections(scat_p.beam_trajectory(1,:), ...
            scat_p.beam_trajectory(3,:), scat_p.beam_trajectory(1,:), p_3dbline, 1);
        if length(p_cross) >= 2
            crossidx = find(p_cross >= scat_p.peak(1), 1);
            if crossidx >= 2
                scat_p.peak_3db = [p_cross(crossidx-1); p_cross(crossidx); ...
                    p_cross(crossidx) - p_cross(crossidx-1)]; % [start, end, angle_range]
            end
        end
    end
    scatterer_points{p} = scat_p;
end
end