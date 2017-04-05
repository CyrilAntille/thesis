%GETIAAMULTIBEAM Multibeam implementation of the iterative adaptive filter (IAA)
%                for linear arrays.
%
% [imSbp imMbp] = getIAAMultiBeam(dataCube, indsI, regCoef, V, beamThetas, scanGridThetas, pitchInLambdas, useMBCovInIterations, nIterations, doCalcSbp, verbose, useSlowForLoop)
%
% dataCube     : A range x beams x array-elements image cube
% indsI        : Indices of ranges that are to be processed; indsI==0 means all ranges
% regCoef      : Regularization coefficient for diagonal loading
% V            : The columns span an orthogonal subspace; if V is empty ([]) then the optimal subspace is estimated from eigenvalue decomposition
% beamThetas     : The angles for the actual beams
% scanGridThetas : The scanning grid in angles (could sample much denser than that of physical beams)
% pitchInLambdas : Distance between array elements in center wavelengths
% useMBCovInIterations: Whether the power estimates in the loop should be based on the closest "physical" beam or the multi beam covariance estimate
% nIterations  : Number of iterations; 5 to 15 iterations seem often adequate
% doCalcSbp : Whether or not to calculate the final power estimates from single-beam values (linearly interpolated at the scanGridThetas)
% useSlowForLoop : {Default: false} Enforces a slow for-loop over the scangrid angles [Perhaps useful for comparisons with approaches not that easitly matlab-optimized?]
% 
%
% Note: Matrix inversion is done like this: Ri = pinv(Rline + regCoef/K*trace(Rline)*I);
%
% Last modified:
% 2009.11.23 - Are C. Jensen {Created the function}
% 2009.12.09 - Are C. Jensen {Changed amplitude output to use linearly interpolated beam responses (instead of nearest neighbor)}
% 2010.06.30 - Are C. Jensen {Made algorithm skip interpolation when beamThetas==scanGridThetas}
% 2011.01.21 - Are C. Jensen {Fixed bug when ~useMultibeamPowerEstimate, and introduced simple error-check for scanGridThetas values  when needs to interpolate.}
% 2011.06.23 - Are C. Jensen {Renamed parameters signalling the use of 'sbp' and 'mpb' instead of the misleading power and amplitude terms, and introduced the useSlowForLoop parameter}
% 2011.11.25 - Are C. Jensen {Switched to beamThetas and scanGridThetas in angles, and added the "pitchInLambdas" parameter }
% 2017.04.05 - Cyril Antille {If V is emtpy, the optimal subspace value is estimated from eigenvalue decomposition}
function [imSbp, imMbp] = getIAAMultiBeam(dataCube, indsI, regCoef, V, beamThetas, scanGridThetas, pitchInLambdas, useMBCovInIterations, nIterations, doCalcSbp, verbose, useSlowForLoop)

if nargin<12
  useSlowForLoop = false;
end

[N, M, K] = size(dataCube);

if indsI<=0, indsI = 1:N; end;
nBeams = M;
doInterpolate = length(scanGridThetas)~=length(beamThetas) || sum(abs(scanGridThetas-beamThetas))>1e-10;

if doCalcSbp && doInterpolate
  if min(scanGridThetas)<min(beamThetas) || max(scanGridThetas)>max(beamThetas)
    warning(' Elements of scanGridThetas should/must not be outside the range of beamThetas.');
  end
end

% Create a matrix, A, having columns representing the scanning grid
scanGridPsis = 2*pi*sin(scanGridThetas)*pitchInLambdas;
nSources = length(scanGridThetas);
A = zeros(K,nSources);
n = (0:K-1)';
for m=1:nSources
  A(:,m) = exp( complex(0,1)*scanGridPsis(m)*n );
end

if isempty(V)
    %V = eye(K); % Fullspace span
    eigenvals = svd(A);
    eigenthld = (max(eigenvals) - min(eigenvals)) / 2;
    nSubDims = find(eigenvals <= eigenthld, 1);
    if isempty(nSubDims)
        nSubDims = K;
    end
    fprintf('IAAMultiBeam nSubDims from eigenvalues: %d\n', nSubDims)
    %V = diag(eigenvals(1:nSubDims));
    V = getSimpleBeamspace(K, nSubDims);
else
    nSubDims = length(V(1,:));
end
I = eye(nSubDims);
indsOneToNSubdims = (1:nSubDims)';
A = V'*A;

% Build the steering matrix, Ab, in which each column represents one "physical" beam direction
beamPsis = 2*pi*sin(beamThetas)*pitchInLambdas;
Ab = zeros(K,nBeams);
n = (0:K-1)';
for b=1:nBeams
  Ab(:,b) = exp( complex(0,1)*beamPsis(b)*n );
end

imSbp = zeros(N,nSources, nIterations);
imMbp = zeros(N,nSources, nIterations);

% For every time instance (depth) in the image, calculate the output power
% estimate for every scan-grid point
for i=indsI
  % "Shift" the beams so as if they were arriving from their original angle
  shiftedArrayResponses = V' * ( Ab.*transpose(squeeze(dataCube(i,:,:))) );

  R = shiftedArrayResponses*shiftedArrayResponses' / nBeams; % 20120221: FB-averaging converges to the same solution..

  if doCalcSbp && doInterpolate
    % Find single-beam array responses at scanGridThetas using linear
    % interpolation
    shiftedArrayResponses = interp2(beamThetas, indsOneToNSubdims, shiftedArrayResponses, scanGridThetas, indsOneToNSubdims, 'linear');
  end

  p = abs(diag( A'*R*A )) / K^2;  % A first estimate of the signal powers, and temp storage of mpb estimates
  imMbp(i,:,1) = p;
  sbpTmp = zeros(nSources,1);   % Temp storage for single-beam-power estimates

  for iter=2:nIterations
    Rline = A*diag(p)*A';  % Estimate of covariance matrix using the sparse model
    Ri = pinv(Rline + regCoef/K*trace(Rline)*I);

    % Whether or not we're forced to use matlab-slow for-loop [Remove this "feature"?]
    if useSlowForLoop
      
      for k=1:nSources
        ak = A(:,k);
        wk = Ri*ak / ( ak'*Ri*ak);
        if doCalcSbp
          sbpTmp(k) = abs(wk'*shiftedArrayResponses(:, k))^2;
        end
        if useMBCovInIterations
          p(k) = abs(wk'*R*wk);
        else
          p(k) = sbpTmp(k);
        end
      end
      imMbp(i,:,iter) = p;
      imSbp(i,:,iter) = sbpTmp;
    
    else
      
      % Find the weights for every scangrid-angle
      W = Ri*A * diag( 1./diag(A'*Ri*A) );
      
      % Calculate the multi-beam-power (MBP) and, if so chosen, the single-beam-power (SBP) outputs
      p = abs(diag(W'*R*W));
      
      imMbp(i,:,iter) = p;
      if doCalcSbp || ~useMBCovInIterations
        imSbp(i,:,iter) = abs(diag(W'*shiftedArrayResponses)).^2';
      end
      
%       if max( abs(p' - imMbp(i,:,iter-1)) ) < 1e-3*max(imMbp(i,:,1))
%         imMbp(i,:,end) = p;
%         imSbp(i,:,end) = p;
%         break;
%       end
      
      if ~useMBCovInIterations
        p = imSbp(i,:,iter);
      end
      
    end
  end

  if verbose
    percent = round( 100*(i-min(indsI))/(max(indsI)-min(indsI)) );
    if mod(i,5)==0, fprintf(' %d%%', percent); end;
    if mod(i,100)==0, fprintf('\n'); end;
  end
end

