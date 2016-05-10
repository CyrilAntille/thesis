%GETCAPONMULTIBEAM Multibeam implementation of the Minimum Variance / Capon
%                  filter for linear arrays.
%
% [imSbp imMbp] = getCaponMultiBeam(dataCube, indsI, regCoef, V, beamPsis, scanGridPsis, returnSbpAndOrMbp, doForwardBackwardAveraging, verbose, useSlowForLoop)
%
% dataCube     : A range x beams x array-elements image cube
% indsI        : Indices of ranges that are to be processed; indsI==0 means all ranges
% regCoef      : Regularization coefficient for diagonal loading
% V            : The columns span an orthogonal subspace; if V is empty ([]) than the full space is used (no subspacing)
% beamPsis     : The angles for the actual beams given in "psi" coordinates ( e.g. beamPsis = pi*sin(beamThetas); )
% scanGridPsis : The scanning grid given in "psi"-space coordinates (could sample much denser than that of physical beams)
% returnSbpAndOrMbp : (1) SBP: Only calculate power estimates from single-beam values (linearly interpolated at the scanGridPsis) 
%                     (2) MBP: Calculate output power based on all beams, i.e., w'*R*w
%                     (3) Do both
% doForwardBackwardAveraging : Whether to do forward-backward averaging
% useSlowForLoop : {Default: false} Enforces a slow for-loop over the scangrid angles [Perhaps useful for comparisons with approaches not that easitly matlab-optimized?]
% 
%
% Note: Matrix inversion is done like this: Ri = inv(R + regCoef/K*trace(R)*I); ('K': number of array elements)
%
% Last modified:
% 2009.11.23 - Are C. Jensen {Created the function}
% 2009.12.09 - Are C. Jensen {Changed amplitude output to use linearly interpolated beam responses (instead of nearest neighbor)}
% 2010.06.30 - Are C. Jensen {Made algorithm skip interpolation when beamPsis==scanGridPsis}
% 2010.07.05 - Are C. Jensen {Added support for forward-backward-averaging}
% 2011.01.21 - Are C. Jensen {Removed parameter doCalcAmplitude and introduced returnPowAndOrAmp}
% 2011.06.23 - Are C. Jensen {Renamed parameters signalling the use of 'sbp' and 'mpb' instead of the misleading power and amplitude terms, and introduced the useSlowForLoop parameter}
function [imSbp imMbp] = getCaponMultiBeam(dataCube, indsI, regCoef, V, beamPsis, scanGridPsis, returnSbpAndOrMbp, doForwardBackwardAveraging, verbose, useSlowForLoop)

if nargin<10
  useSlowForLoop = false;
end

[N M K] = size(dataCube);

if indsI<=0, indsI = 1:N; end;
if isempty(V), V = eye(K); end;
nBeams = M;
nSubDims = length(V(1,:));
I = eye(nSubDims);
doInterpolate = length(scanGridPsis)~=length(beamPsis) || sum(abs(scanGridPsis-beamPsis))>1e-10;
indsOneToNSubdims = (1:nSubDims)';

doCalcSbp = false;
doCalcMbp = false;
if returnSbpAndOrMbp==1 || returnSbpAndOrMbp==3
  doCalcSbp = true;
end
if returnSbpAndOrMbp==2 || returnSbpAndOrMbp==3
  doCalcMbp = true;
end

% Create a matrix, A, having columns representing the scanning grid
nSources = length(scanGridPsis);
A = zeros(K,nSources);
n = (0:K-1)';
for m=1:nSources
  A(:,m) = exp( complex(0,1)*scanGridPsis(m)*n );
end
A = V'*A;

% Build the steering matrix, Ab, in which each column represents one "real-life" beam direction
Ab = zeros(K,nBeams);
n = (0:K-1)';
for b=1:nBeams
  Ab(:,b) = exp( complex(0,1)*beamPsis(b)*n );
end

imSbp = zeros(N,nSources);
imMbp = zeros(N,nSources);

% For every time instance (depth) in the image, calculate the output power
% estimate for every scan-grid point
for i=indsI
  % "Shift" the beams so as if they were all directed perpendicular to the array
  shiftedArrayResponses = ( Ab.*transpose(squeeze(dataCube(i,:,:))) );

  if doForwardBackwardAveraging
    % If we _are_ using FB averaging, we have to wait until after forming
    % the covariance matrix to go to the reduced-size space:
    R = shiftedArrayResponses*shiftedArrayResponses' / nBeams;
    R = 0.5*(R + rot90(conj(R), 2));
    R = V'*R*V;
    if doCalcSbp
      shiftedArrayResponses = V'*shiftedArrayResponses;
    end
  else
    shiftedArrayResponses = V'*shiftedArrayResponses;
    R = shiftedArrayResponses*shiftedArrayResponses' / nBeams;
  end

  if doCalcSbp && doInterpolate
    % Find single-beam array responses at scanGridPsis using linear
    % interpolation
    shiftedArrayResponses = interp2(beamPsis, indsOneToNSubdims, shiftedArrayResponses, scanGridPsis, indsOneToNSubdims, 'linear');
  end

  Ri = inv(R + regCoef/K*trace(R)*I);
  
  % Whether or not we're forced to use matlab-slow for-loop [Remove this "feature"?]
  if useSlowForLoop
    
    % Loop through the scan-grid points
    for k=1:nSources
      ak = A(:,k);
      w = Ri*ak / (ak'*Ri*ak);
      if doCalcMbp
        imMbp(i,k) = abs(w'*R*w);
      end
      if doCalcSbp
        imSbp(i,k) = abs(w'*shiftedArrayResponses(:, k)).^2;  % Will take abs and square these values at the end of the script
      end
    end
    
  else
    
    % Find the weights for every scangrid-angle
    W = Ri*A * diag( 1./diag(A'*Ri*A) );
    
    % Calculate the multi-beam-power (MBP) and/or the single-beam-power (SBP) outputs
    if doCalcMbp
      imMbp(i,:) = abs(diag(W'*R*W));
    end
    if doCalcSbp
      imSbp(i,:) = abs(diag(W'*shiftedArrayResponses)).^2';
    end
    
  end

  if verbose
    percent = round( 100*(i-min(indsI))/(max(indsI)-min(indsI)) );
    if mod(i,5)==0, fprintf(' %d%%', percent);, end;
    if mod(i,100)==0, fprintf('\n');, end;
  end
end

