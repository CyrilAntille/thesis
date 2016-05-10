% You need to download FieldII and add path to correct directory
%addpath '/hom/dsb/projects/matlab/field/' -end
%addpath '\\tsclient\C\Users\Cyril\Documents\MATLAB\Field_II_ver_3_24' -end
addpath 'C:\Users\Cyril\Documents\MATLAB\Field_II_ver_3_24' -end
  
  field_init(0);
  P = Parameters();
  %points = [40*sind(-1) 0 40*cosd(-1); 0 0 40; 40*sind(2) 0 40*cosd(2)]*1e-3;
  %points = [sind(-1) 0 cosd(-1); sind(2) 0 cosd(2)]*40e-3;
  points = [sind(-0.5) 0 cosd(-0.5); sind(2.5) 0 cosd(2.5)]*40e-3;
  
  Phantom = PointPhantom(points);

  Data = CalcRespAll(P,Phantom);
  field_end();

  save -v7.3 2.5PointsData  points P Data;

  % Beamform data
  BFData = BeamformAll(P,Data);

  save -v7.3 2.5PointsBeamFormed  P BFData;

