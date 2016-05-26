classdef Shift
    %SHIFT Defines shift per frame
    properties
        type = ShiftType.RadialVar;
        val = 0; % Ref ShiftType
        num_shifts = 0;
    end
    methods
        function obj = Shift(type, val, num_shift)
            obj.type = type;
            obj.val = val;
            obj.num_shifts = num_shift;
        end
        function shifts = getShifts(obj, P)
            shift_val = obj.val;
            if obj.type == ShiftType.RadialVar
                shift_val = (P.Tx.Theta(2) - P.Tx.Theta(1)) * shift_val;
            elseif obj.type == ShiftType.LateralVar
                shift_val = (P.Tx.Theta(2) - P.Tx.Theta(1)) * shift_val;
                shift_val = sin(shift_val) * P.Tx.FocRad;
            elseif obj.type == ShiftType.LateralSpeed
                shift_val = obj.val * 1e-4; % 0.1ms per beam
            elseif obj.type == ShiftType.LateralSpeedFrame
                shift_val = obj.val * 1e-4 * P.Tx.NTheta; % 0.1ms per beam
            end
            shifts = (0:obj.num_shifts-1) * shift_val;
        end
        function s_phantom = shiftPositions(obj, phantom, shift_val)
            s_phantom = copyStruct(phantom);
            if obj.type == ShiftType.RadialVar || ...
                    obj.type == ShiftType.RadialCst
                [Theta,Phi,R] = cart2sph(s_phantom.positions(:,3),...
                    s_phantom.positions(:,1), s_phantom.positions(:,2));
                Theta = Theta + shift_val;
                [z,x,y] = sph2cart(Theta,Phi,R);
                s_phantom.positions = [x,y,z];
            else %linear shift (along x)
                s_phantom.positions = [s_phantom.positions(:,1) + shift_val,...
                    s_phantom.positions(:,2), s_phantom.positions(:,3)];
            end
        end
    end
end

