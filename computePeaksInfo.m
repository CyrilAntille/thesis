function [ data_peaks ] = computePeaksInfo( mainP, data_phantom, ...
    data_DA, data_BF )
%COMPUTEPEAKSINFO Summary of this function goes here
%   Detailed explanation goes here

data_peaks = cell([1, length(mainP.methods_set)]);
for m=1:length(mainP.methods_set)
    m_BF = data_BF{m};
    m_P = mainP.P;
    m_shift = mainP.shift;
    if strcmp(mainP.methods_set(m),'IAA-MBMB-Upsampled')
        m_P = mainP.copyP(mainP.upsample_number);
        ratio = mainP.num_beams / m_P.Tx.NTheta;
        m_shift = Shift(mainP.shift.type, ratio * mainP.shift.val, ...
            m_P.Tx.NTheta, mainP.shift.direction);
    end
    shifts = m_shift.getShifts(m_P);
    
    bf_img = db(m_BF);
    scatterer_points = cell([1 length(mainP.pts_range)]);
    pts_center = zeros([1 length(mainP.pts_range)]);
    % Note: trajectory same for all BFs except IAA-Upsampled
    for p=1:length(mainP.pts_range)
        scat_p = struct;
        scat_p.p_trajectory = zeros([3 m_shift.num_shifts]); % [x, z, ampl]
        scat_p.beam_trajectory = zeros([3 m_shift.num_shifts]);
        for s=1:m_shift.num_shifts
            s_phantom = m_shift.shiftPositions(data_phantom, shifts(s));
            scat_p.p_trajectory(:,s) = [s_phantom.positions(p, 1), ...
                s_phantom.positions(p, 3), s_phantom.amplitudes(p)];

            s_az = tan(m_P.Tx.Theta(s)) * s_phantom.positions(p, 3);
            s_radius = sqrt(s_phantom.positions(p, 3)^2 + s_az^2);
            s_radius_idx = find(data_DA.Radius >= s_radius,1);
            if ~isempty(s_radius_idx)
                scat_p.beam_trajectory(:,s) = [s_az; s_radius; ...
                    bf_img(s_radius_idx, s)];
            end
        end
        % Expected point azimuth when beam hits it
        diff_az = abs(scat_p.p_trajectory(1,:) - scat_p.beam_trajectory(1,:));
        [~, az_idx] = min(diff_az);
        scat_p.pt_center = scat_p.beam_trajectory(1, az_idx);
        pts_center(p) = scat_p.pt_center;
        scatterer_points{p} = scat_p;
    end
    
    for p=1:length(mainP.pts_range)
        scat_p = scatterer_points{p};
        % Finds nearest peak
        [pdb, paz] = findpeaks(scat_p.beam_trajectory(3,:), ...
            scat_p.beam_trajectory(1,:));

        [~, paz_idx] = min(abs(paz - scat_p.pt_center));
        [~, pt_idx] = min(abs(pts_center - paz(paz_idx)));
        if scat_p.pt_center - paz(paz_idx) < 1e-3 && pt_idx == p
            if paz(paz_idx) - scat_p.pt_center < 0 && paz_idx < length(paz)
                % -> compare with next index
                [~, pt_idx] = min(abs(pts_center - paz(paz_idx+1)));
                if pdb(paz_idx+1) > pdb(paz_idx) && pt_idx == p && ...
                        pt_center - paz(paz_idx + 1) < 1e-3
                    paz_idx = paz_idx + 1;
                end
            elseif paz(paz_idx) - scat_p.pt_center >= 0 && paz_idx > 0
                % -> compare with previous index
                [~, pt_idx] = min(abs(pts_center - paz(paz_idx-1)));
                if pdb(paz_idx-1) > pdb(paz_idx) && pt_idx == p && ...
                        pt_center - paz(paz_idx - 1) < 1e-3
                    paz_idx = paz_idx - 1;
                end
            end
            scat_p.peak = [paz(paz_idx); pdb(paz_idx)];
            % Peak 3dB width
            p_3dbline = ones([1 size(scat_p.p_trajectory,2)]) .* (scat_p.peak(2) - 3);
            [p_x, ~, ~, ~] = intersections(scat_p.beam_trajectory(1,:), ...
                scat_p.beam_trajectory(3,:), scat_p.beam_trajectory(1,:), p_3dbline, 1);
            if length(p_x) >= 2
                p_x2 = find(p_x >= scat_p.peak(1), 1);
                if p_x2 >= 2
                    scat_p.peak_3db = [p_x(p_x2-1); p_x(p_x2); ...
                        p_x(p_x2) - p_x(p_x2-1)];
                end
            end
        end
        scatterer_points{p} = scat_p;
    end
    data_peaks{m} = scatterer_points;
end

