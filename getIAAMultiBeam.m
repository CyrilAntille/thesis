%GETIAAMULTIBEAM Multibeam implementation of the iterative adaptive filter (IAA)
%                for linear arrays.
%
% [imAmplitude imPower] = getIAAMultiBeam(dataCube, indsI, regCoef, V, beamPsis, scanGridPsis, useMultibeamPowerEstimate, nIterations, doCalcAmplitude, verbose)
%
% dataCube     : A range x beams x array-elements image cube
% indsI        : Indices of ranges that are to be processed; indsI==0 means all ranges
% regCoef      : Regularization coefficient for diagonal loading
% V            : The columns span an orthogonal subspace; if V is empty ([]) than the full space is used (no subspacing)
% beamPsis     : The angles for the actual beams given in "psi" coordinates ( e.g. beamPsis = pi*sin(beamThetas); )
% scanGridPsis : The scanning grid given in "psi"-space coordinates (could sample much denser than that of physical beams)
% useMultibeamPowerEstimate : Whether the power estimates in the loop should be based on the closest "physical" beam or the multi beam covariance estimate
% nIterations  : Number of iterations; 5 to 15 iterations seem often adequate
% doCalcAmplitude : Whether or not to calculate amplitude estimates from single-beam values (linearly interpolated at the scanGridPsis)
% 
%
% Note: Matrix inversion is done like this: Ri = pinv(Rline + regCoef/K*trace(Rline)*I);
%
% Last modified:
% 2009.11.23 - Are C. Jensen {Created the function}
% 2009.12.09 - Are C. Jensen {Changed amplitude output to use linearly interpolated beam responses (instead of nearest neighbor)}
% 2010.06.30 - Are C. Jensen {Made algorithm skip interpolation when beamPsis==scanGridPsis}
function [imAmplitude imPower] = getIAAMultiBeam(dataCube, indsI, regCoef, V, beamPsis, scanGridPsis, useMultibeamPowerEstimate, nIterations, doCalcAmplitude, verbose)

[N M K] = size(dataCube);

if indsI<=0, indsI = 1:N;, end;
if isempty(V), V = eye(K);, end;
nBeams = M;
nSubDims = length(V(1,:));
I = eye(nSubDims);
doInterpolate = length(scanGridPsis)~=length(beamPsis) || sum(abs(scanGridPsis-beamPsis))>1e-15;

% Create a matrix, A, having columns representing the scanning grid
nSources = length(scanGridPsis);
A = zeros(K,nSources);
n = (0:K-1)';
for m=1:nSources
  A(:,m) = exp( complex(0,1)*scanGridPsis(m)*n );
end
A = V'*A;

% Build the steering matrix, Ab, in which each column represents one "physical" beam direction
Ab = zeros(K,nBeams);
n = (0:K-1)';
for b=1:nBeams
  Ab(:,b) = exp( complex(0,1)*beamPsis(b)*n );
end

imAmplitude = zeros(N,nSources, nIterations);
imPower = zeros(N,nSources, nIterations);

% For every time instance (depth) in the image, calculate the output power
% and amplitude estimate for every scan-grid point
for i=indsI
  % "Shift" the beams so as if they were all directed perpendicular to the array
  shiftedArrayResponses = V' * ( Ab.*transpose(squeeze(dataCube(i,:,:))) );

  if doCalcAmplitude && doInterpolate
    % Find single-beam array responses at scanGridPsis using linear
    % interpolation
    shiftedArrayResponsesInterpolatedScangrid = interp2(beamPsis, (1:nSubDims)', shiftedArrayResponses, scanGridPsis, (1:nSubDims)', 'linear');
  end

  R = shiftedArrayResponses*shiftedArrayResponses' / nBeams;

  p = diag( A'*R*A ) / K^2;  % A first estimate of the signal powers
  imPower(i,:,1) = p;
  amp = zeros(nSources,1);   % Temp storage for single-beam-amplitude estimates

  for iter=2:nIterations
    Rline = A*diag(p)*A';  % Estimate of covariance matrix using the sparse model
    Ri = pinv(Rline + regCoef/K*trace(Rline)*I);
     
    for k=1:nSources
      ak = A(:,k);
      wk = Ri*ak / ( ak'*Ri*ak);
      if doCalcAmplitude
        if doInterpolate
          amp(k) = wk'*shiftedArrayResponsesInterpolatedScangrid(:, k);
        else
          amp(k) = wk'*shiftedArrayResponses(:, k);
        end
      end
      if useMultibeamPowerEstimate
        p(k) = abs(wk'*R*wk);
      else
        p(k) = abs(imAmplitude(i,k,iter)).^2;
      end
    end
    imPower(i,:,iter) = p;
    imAmplitude(i,:,iter) = amp;
  end

  if verbose
    percent = round( 100*(i-min(indsI))/(max(indsI)-min(indsI)) );
    if mod(i,5)==0, fprintf(' %d%%', percent);, end;
    if mod(i,100)==0, fprintf('\n');, end;
  end
end

