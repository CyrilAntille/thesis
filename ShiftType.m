classdef ShiftType < handle
    %SHIFTTYPE Defines shift per frame
    properties
        shift = 0
        num_shifts = 0
    end
    methods
        function shifts = getShifts(obj, P)
            val = obj.shift;
            if obj == ShiftType.RadialVar
                val = (P.Tx.Theta(2) - P.Tx.Theta(1)) * val;
            elseif obj == ShiftType.LateralVar
                val = (P.Tx.Theta(2) - P.Tx.Theta(1)) * val;
                val = sin(val) * P.Tx.FocRad;
            end
            shifts = (0:obj.num_shifts-1) * val;
        end
        function s_phantom = shiftPositions(obj, phantom, shift_val)
            s_phantom = copyStruct(phantom);
            if obj == ShiftType.RadialVar || obj == ShiftType.RadialCst
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
    enumeration
        LateralCst %in meters
        RadialCst % in radians
        RadialVar % in beam distance ratio(1 = moves from one beam to another).
        LateralVar % Same. Lateral shift calculated at focus range.
    end
end