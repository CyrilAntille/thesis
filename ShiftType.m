classdef ShiftType
    %SHIFTTYPE Defines shift type
    enumeration
        LinearCst %in meters per shift
        RadialCst % in radians per shift
        RadialVar % in beam distance ratio per shift
        % (1 = moves from one beam to another).
        LinearVar % Same. Shift value calculated at focus range.
        
        % Assumes max range of 15cm (-> 30cm two-way) -> beam transmit
        % lasts for 300mm / c (1500 m/s) = 2 * 1e-4 seconds = 0.2 ms
        LinearSpeed % in meters/seconds
        LinearSpeedFrame % Same as LateralSpeed multiplied by num beams
    end
    methods(Static)
        function unit = getShiftTypeUnit(type)
            unit = 'beam distance ratio / shift';
            if type == ShiftType.LinearSpeed || ...
                    type == ShiftType.LinearSpeedFrame
                unit = 'meter/second';
            elseif type == ShiftType.LinearCst
                unit = 'meter/shift';
            elseif type == ShiftType.RadialCst
                unit = 'radian/shift';
            end
        end
    end
end