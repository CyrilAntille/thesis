classdef Shift
    %SHIFT Defines shift per frame/beam
    properties
        type = ShiftType.RadialVar;
        val = 0; % Ref ShiftType.m
        num_shifts = 0;
        direction = 0; % Motion direction, in degrees.
        % Values must be between -180 and 180 degrees. 0 = left to right
        % motion, 90 = moving away from array, 180 (=-180) = right to left.
        stop_and_go_ratio = 1; % 'number of beams' per shift
        % If stop_and_go_ratio = 1, then the motion is continuous.
        % If > 1, then stop-and-go. Example: = 3 -> shift once every 3
        % beams. Note that the shift is of the same value as if continuous.
    end
    methods
        function obj = Shift(type, val, num_shift, direction, stop_and_go_ratio)
            if nargin == 4
                stop_and_go_ratio = 1;
            end
            obj.type = type;
            obj.val = val;
            obj.num_shifts = num_shift;
            obj.direction = direction;
            obj.stop_and_go_ratio = stop_and_go_ratio;
        end
        function shifts = getShifts(obj, P)
            shift_val = obj.val;
            if obj.type == ShiftType.RadialVar
%                 shift_val = (P.Tx.Theta(2) - P.Tx.Theta(1)) * shift_val;
                shift_val = (P.Tx.SinTheta(2) - P.Tx.SinTheta(1)) * shift_val;
            elseif obj.type == ShiftType.LinearVar
                shift_val = (P.Tx.SinTheta(2) - P.Tx.SinTheta(1)) * ...
                    shift_val * P.Tx.FocRad;
            elseif obj.type == ShiftType.LinearSpeed
                shift_val = obj.val * 2*1e-4; % 0.2ms per beam
            elseif obj.type == ShiftType.LinearSpeedFrame
                shift_val = obj.val * 2*1e-4 * P.Tx.NTheta; % 0.2ms per beam
            end
            max_shift = floor(obj.num_shifts/(2*obj.stop_and_go_ratio));
            shifts = (-max_shift:max_shift) * shift_val;
            if mod(obj.num_shifts, 2) == 0
                shifts = shifts(2:end);
            end
            if obj.type == ShiftType.RadialVar || obj.type == ShiftType.LinearVar
                shifts = asin(shifts);
            end
            shifts = repelem(shifts, obj.stop_and_go_ratio); % stop-and-go
        end
        function s_phantom = shiftPositions(obj, phantom, shift_val)
            % Note: - shift_val = shift from right to left (+ from left...)
            s_phantom = copyStruct(phantom);
            if obj.type == ShiftType.RadialVar || ...
                    obj.type == ShiftType.RadialCst
                [Theta,Phi,R] = cart2sph(s_phantom.positions(:,3),...
                    s_phantom.positions(:,1), s_phantom.positions(:,2));
                if obj.direction > -90 && obj.direction < 90
                    % For radial motion, the only two possible directions
                    % are left to right (0) or right to left (-180).
                    Theta = Theta + shift_val;
                else
                    Theta = Theta - shift_val;
                end
                [z,x,y] = sph2cart(Theta,Phi,R);
                s_phantom.positions = [x,y,z];
                % New approach
%                 for p=1:size(s_phantom.positions, 1)
%                     pos = num2cell(s_phantom.positions(p, :));
%                     [x, y, z] = deal(pos{:});
%                     Theta = atan(x/z);
%                     R = z / cos(Theta);
%                     if obj.direction > -90 && obj.direction < 90
%                     % For radial motion, the only two possible directions
%                     % are left to right (0) or right to left (-180).
%                         Theta = Theta + shift_val;
%                     else
%                         Theta = Theta - shift_val;
%                     end
%                     x = sin(Theta) * R;
%                     z = cos(Theta) * R;
%                     s_phantom.positions(p,:) = [x,y,z];
%                 end
            else %linear shift (along obj.direction)
                s_phantom.positions = [s_phantom.positions(:,1) ...
                    + shift_val * cosd(obj.direction), ...
                    s_phantom.positions(:,2), s_phantom.positions(:,3) ...
                    + shift_val * sind(obj.direction)];
            end
        end
    end
end

