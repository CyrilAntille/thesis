% -----------------------------------------------------------------
%
% Function     : function Phantom = PointPhantom();
%
% Description  : Defines a simple phantom.
%
%%

function Phantom = PointPhantom(points, ampl)
% Define fantom
  % input:  points : matrix with phantom postitons [x y z]
  
Phantom.positions = points;
Phantom.amplitudes = ones(size(points, 1)) * ampl;
  
