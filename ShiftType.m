classdef ShiftType < handle
    %SHIFTTYPE Defines shift type
    enumeration
        LateralCst %in meters
        RadialCst % in radians
        RadialVar % in beam distance ratio(1 = moves from one beam to another).
        LateralVar % Same. Lateral shift calculated at focus range.
    end
end