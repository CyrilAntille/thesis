classdef ShiftType
    %SHIFTTYPE Defines shift type
    enumeration
        LateralCst %in meters
        RadialCst % in radians
        RadialVar % in beam distance ratio(1 = moves from one beam to another).
        LateralVar % Same. Lateral shift calculated at focus range.
        
        % Assumes max range of 15cm (-> 30cm two-way) -> beam transmit
        % lasts for 300mm / c (1500 m/s) = 2 * 1e-4 seconds = 0.2 ms
        LateralSpeed % in meters/seconds
        LateralSpeedFrame % Same as LateralSpeed multiplied by num beams
    end
end