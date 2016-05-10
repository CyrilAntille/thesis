% -----------------------------------------------------------------
%
% Function     : function Phantom = SimplePhantom();
%
% Description  : Defines a simple phantom.
%
%%

function Phantom = SimplePhantom(points, ampl);
% Define fantom
  % input:  points : matrix with phantom postitons [x y z]
  
Phantom.positions = points;
Phantom.amplitudes = ones(size(points, ampl));
  
