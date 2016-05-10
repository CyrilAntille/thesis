function [ bf_im ] = ComputeBF( BFData, P, bf_method )
%COMPUTEBF Summary of this function goes here
%   Detailed explanation goes here
dl = 5/100; % Diagonal loading
sbl = P.Rx.no_elements/2; % Subarray length (48)
dofb = true; % Do forward-backward averaging
% nsd = ceil(rad2deg(P.Tx.Theta(end))); % Number of subdimensions (to compute V)
nsd = 25; % Number of subdimensions (to compute V)
dataCube = permute(BFData,[2,1,3]);
V = getSimpleBeamspace(P.Rx.no_elements,nsd);
if bf_method == 1 %  MV
    warning('off');
    MVImage = getCapon(dataCube,0,0,dl,sbl,2,5,0,dofb,0);
    warning('on');
    bf_im = MVImage;
%     bf_im = 20*log10(abs(MVImage));
elseif bf_method == 2 % MV-Multibeam
    MVMultiImage = getCaponMultiBeam(dataCube,0,dl,V, ...
        pi*P.Tx.SinTheta,pi*P.Tx.SinTheta,1,0,1,0);
    bf_im = MVMultiImage;
%     bf_im = 10*log10(abs(MVMultiImage));
elseif bf_method == 3 % IAA
    [IAAImageAmp, IAAImagePow] = getIAAMultiBeam(dataCube,0,dl,V, ...
        pi*P.Tx.SinTheta,pi*P.Tx.SinTheta,1,10,1,1);
    bf_im = IAAImageAmp(:,:,10);
%     bf_im = 20*log10(abs(IAAImageAmp(:,:,10)));
elseif bf_method == 4 % IAA 500 pts
    % NOTE: Should scanGridPsis be uniform along sin(theta) or theta???
    NTheta = 500;
    SinThMax = 0.3; 
    scanGridPsis = pi*linspace(-SinThMax,SinThMax,NTheta);
    [IAAImageAmp, IAAImagePow] = getIAAMultiBeam(dataCube,0,dl,V, ...
        pi*P.Tx.SinTheta,scanGridPsis,1,10,1,1);
    bf_im = IAAImageAmp(:,:,10);
%     bf_im = 20*log10(abs(IAAImageAmp(:,:,10)));
else % DAS
    bf_im = mean(BFData,3);
%     bf_im = 20*log10(abs(mean(BFData,3))).'; % sum of weights equal to one'
end

end

