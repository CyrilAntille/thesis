classdef Shift
    %SHIFT Defines shift per frame/beam
    properties
        type = ShiftType.RadialVar;
        val = 0; % Ref ShiftType.m
        num_shifts = 0;
        direction = 0; % Motion direction, in degrees.
        % Values must be between -180 and 180 degrees. 0 = left to right
        % motion, 90 = moving away from array, 180 (=-180) = right to left.
    end
    methods
        function obj = Shift(type, val, num_shift, direction)
            obj.type = type;
            obj.val = val;
            obj.num_shifts = num_shift;
            obj.direction = direction;
        end
        function shifts = getShifts(obj, P)
            shift_val = obj.val;
            if obj.type == ShiftType.RadialVar
                shift_val = (P.Tx.Theta(2) - P.Tx.Theta(1)) * shift_val;
            elseif obj.type == ShiftType.LinearVar
                shift_val = (P.Tx.Theta(2) - P.Tx.Theta(1)) * shift_val;
                shift_val = sin(shift_val) * P.Tx.FocRad;
            elseif obj.type == ShiftType.LinearSpeed
                shift_val = obj.val * 2*1e-4; % 0.2ms per beam
            elseif obj.type == ShiftType.LinearSpeedFrame
                shift_val = obj.val * 2*1e-4 * P.Tx.NTheta; % 0.2ms per beam
            end
%             shifts = (0:obj.num_shifts-1) * shift_val;
            shifts = (-floor(obj.num_shifts/2):floor(obj.num_shifts/2))...
                * shift_val;
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
            else %linear shift (along x)
                s_phantom.positions = [s_phantom.positions(:,1) ...
                    + shift_val * cosd(obj.direction), ...
                    s_phantom.positions(:,2), s_phantom.positions(:,3) ...
                    + shift_val * sind(obj.direction)];
            end
        end
    end
end

