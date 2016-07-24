function [ bf_im ] = ComputeBF( BFData, P, bf_method )
%COMPUTEBF Summary of this function goes here
%   Detailed explanation goes here
dl = 5/100; % Diagonal loading
sbl = P.Rx.no_elements/2; % Subarray length (48)
dofb = true; % Do forward-backward averaging
% nsd = ceil(rad2deg(P.Tx.Theta(end))); % Number of subdimensions (to compute V)
nsd = 30; % Number of subdimensions (to compute V)
% nsd = 25; % Number of subdimensions (to compute V)
lambda = P.c/P.fc; % For IAA, to get pitchInLambdas
pitchInLambdas = P.Tx.pitch / lambda;
dataCube = permute(BFData,[2,1,3]);
V = getSimpleBeamspace(P.Rx.no_elements,nsd);
if bf_method == 1 %  MV
    warning('off');
    bf_im = getCapon(dataCube,0,0,dl,sbl,2,5,0,dofb,0);
    bf_im = abs(bf_im);
    warning('on');
elseif bf_method == 2 % MV-Multibeam
    bf_im = getCaponMultiBeam(dataCube,0,dl,V, ...
        pi*P.Tx.SinTheta,pi*P.Tx.SinTheta,1,0,1,0);
    bf_im = abs(bf_im);
elseif bf_method == 3 % IAA MB/SB
    [IAAImageAmp, IAAImagePow] = getIAAMultiBeam(dataCube,0,dl,V, ...
        P.Tx.Theta,P.Tx.Theta,pitchInLambdas,1,10,1,1,0);
    bf_im = sqrt(IAAImageAmp(:,:,10));
elseif bf_method == 4 % IAA MB/MB
    [IAAImageAmp, IAAImagePow] = getIAAMultiBeam(dataCube,0,dl,V, ...
        P.Tx.Theta,P.Tx.Theta,pitchInLambdas,1,10,0,1,0);
    bf_im = sqrt(IAAImagePow(:,:,10));
elseif bf_method == 5 % IAA MB/MB oversampled by phase shift
    NTheta = 500;
    scanGridSin = linspace(P.Tx.SinTheta(1),P.Tx.SinTheta(end),NTheta);
    scanGridTheta = asin(scanGridSin);
    [IAAImageAmp, IAAImagePow] = getIAAMultiBeam(dataCube,0,dl,V, ...
        P.Tx.Theta,scanGridTheta,pitchInLambdas,1,10,0,1,0);
    bf_im = sqrt(IAAImagePow(:,:,10));
else % DAS
    bf_im = abs(mean(BFData,3).');
end

end

