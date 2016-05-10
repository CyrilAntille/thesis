%
% $Id:$
%
% -----------------------------------------------------------------
%
% Function     : function Phantom = PlanewaveVesselPhantom(P,VesselAmpl,N,[Seed]);
%
%  OBS: HACKET!!!!!
%
% Description  : Defines the vessel-phantom aka Montaldo, TUFFC, Sept 2009.
%
% Input        : 
%
%   P              : Parameter struct P defined in Parameters.m
%   VesselAmpl     : Amplitude (in dB) of the
%                    vessel areas relative to 0dB.
%   N              : NumberOfPoints in Cyst.
%   Seed           : Seed for initialization of random generator (optional)
%
% Language     : Matlab R2010a
%
% Written by   : Andreas Austeng, Ifi-UiO
%
% Version no.  : 1.0  AA 28.04.11 first version
%
% -----------------------------------------------------------------
% Description:
%
%  OBS: HACKET!!!!!
%
% The cyst is defines using cartesian and spherical coordinates. Points
% are distributed uniformly in cartesian coordinates, but truncated in R,
% Theta and Phi.
%
% The vessels are almost as in Montaldo, TUFFC March 2009, i.e
% One column of vessels with diameter 2mm at x-coord -10mm and z-coords [20
% 30 40 50 60], one column of vessels with d=4mm at x-coord 0mm, same
% z-coords, and one column of vessels with d=8mm at x-coord 10mm, same
% z-coords.
%
% The height of the phantom is defined as 2 times the mainlobe-width
% (covering approximately the main lobe and the first sidelobe on each
% side.
%
% OBS: NO POINTS
% Two pairs of bright points are also defined. The amplitude of these are
% set to 40dB above the speckle level
%
%  OBS: HACKET!!!!!
% -----------------------------------------------------------------

function [Phantom] = PlanewaveVesselPhantom(P,VesselAmpl,N,Seed)
  if ~isdeployed
  help PlanewaveVesselPhantom;
  end

  if nargin<4, Seed=NaN; end
  
% Define cyct-dimensions

Rmin = 5e-3;
Rmax = 75e-3;
ThetaMin = -90 / 180*pi;  % Azimuth angles
ThetaMax = -ThetaMin;
% PhiMin = - asin ( 2*pi/(P.Tx.height/(P.c/P.fc)));
% PhiMax = -PhiMin;
PhiMin = 0;
PhiMax = 0;

% Set seed
if ~isnan(Seed),
    defaultStream = RandStream.getGlobalStream;
    reset(defaultStream,Seed);
end

% Random points
% Define extend of large phantom 
Zmin = 0; Zmax = Rmax;
Xmin = Rmax*sin(ThetaMin); Xmax = Rmax*sin(ThetaMax);
Ymin = Rmax*sin(PhiMin); Ymax = Rmax*sin(PhiMax);
% Ymin = 0; Ymax =0;

%%

Zs = rand(N,1)*(Zmax-Zmin) + Zmin;
Xs = rand(N,1)*(Xmax-Xmin) + Xmin;
Ys = rand(N,1)*(Ymax-Ymin) + Ymin;

% Convert to polar coordinates 
[Theta,Phi,R] = cart2sph(Zs,Xs,Ys);

% and keep points with
% ThetaMin<Theta < ThetaMax, PhiMin < Phi < PhiMax and
% Rmin < R < Rmax

KoordToKeep = (Theta > ThetaMin & Theta < ThetaMax & ...
                    Phi >= PhiMin & Phi<=PhiMax & ...
                    R > Rmin & R < Rmax);
Theta=Theta(KoordToKeep);
Phi=Phi(KoordToKeep);
R=R(KoordToKeep);

% before back-conversion to Cartesian coordinates.

[z,x,y] = sph2cart(Theta,Phi,R);

% Set all amplitudes initially to 1
amps = ones(size(z));

% Set amplitudes in different circular areas:
% The areas lies along a circular path with center -3, -1.5, 1.5 and 3 cm
% off broadside.

zCoords = [20 30 40 50 60]*1e-3; 
xCoords = [-10 0 10]*1e-3;
Diams   = [2 4 8]*1e-3;


for ix = 1:length(xCoords)
    for iz=1:length(zCoords)
        inside = ((x - xCoords(ix)).^2 + ...
            (z - zCoords(iz)).^2 < (Diams(ix)/2)^2);
        amps(inside) = 10^(VesselAmpl/20);
    end
end
% Define one set of double points by moving the two last random points
sep = 3.5e-3; 
Zpoint = 50e-3;

%ii=N-1;
%x(ii) = -sep/2; x(ii+1) = sep/2;
%y(ii) = 0; y(ii+1) = 0;
%z(ii) = Zpoint ; z(ii+1)=z(ii);
%amps(ii) = 10^(40/20); amps(ii+1) = amps(ii);

% Return variables
I = find(amps>0);
Phantom.positions = [x(I) y(I) z(I)];
Phantom.amplitudes = amps(I);
  

