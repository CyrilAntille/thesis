%GETSIMPLEBEAMSPACE Returns a matrix which columns span the
%                   'central-pointing' beamspace.  Used in conjunction with
%                   getApes and getCapon.
%
% [V] = getSimpleBeamspace(L, nSubDimensions)
%
% L              : Dimensionality of the original space
% nSubDimensions : Dimensionality of the subspace
%
%
% Last modified:
% 2009.08.25 - Are C. Jensen {Created the function}
function [V] = getSimpleBeamspace(L, nSubDimensions)

n=0:L-1;
j = complex(0,1);
phi = 0;

V = zeros(L,nSubDimensions);
for m=0:nSubDimensions-1
  V(:,m+1) = exp(-j*(n-(L-1)/2)*(phi-(m-floor(nSubDimensions/2))*2*pi/L));
end

% Normalize the columns of V
for i=1:nSubDimensions
  V(:,i) = V(:,i) / norm(V(:,i));
end
