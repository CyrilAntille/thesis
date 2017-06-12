function [ bf_im ] = ComputeBF( BFData, mainP, bf_method, verbose )
%COMPUTEBF Computes beamformed image with chosen beamformer

lambda = mainP.P.c/mainP.P.fc; % For IAA, to get pitchInLambdas
pitchInLambdas = mainP.P.Tx.pitch / lambda;
dataCube = permute(BFData,[2,1,3]);
V = getSimpleBeamspace(mainP.P.Rx.no_elements,mainP.nsd);
scanGridTheta = mainP.getScanGrid(bf_method); % = mainP.P.Rx.Theta if no upsampling
if strcmp(bf_method, 'DAS')
    bf_im = abs(mean(BFData,3).');
elseif strcmp(bf_method, 'MV')
    warning('off');
    sbl = mainP.sbl * mainP.P.Rx.no_elements;
    bf_im = getCapon(dataCube,0,0,mainP.dl,sbl,mainP.tavg,[],mainP.dofb,verbose);
    bf_im = abs(bf_im);
    warning('on');
elseif strcmp(bf_method, 'MV-MB')
    bf_im = getCaponMultiBeam(dataCube,0,mainP.dl,V, ...
        pi*mainP.P.Rx.SinTheta,pi*mainP.P.Rx.SinTheta,1,0,verbose,0);
    bf_im = abs(bf_im);
elseif strfind(bf_method, 'MV-')
    sbl = strsplit(bf_method, '-');
    sbl = round(str2double(sbl(end)) * mainP.P.Rx.no_elements);
    if sbl >= mainP.P.Rx.no_elements
        sbl = mainP.P.Rx.no_elements - 1;
    end
    warning('off');
    bf_im = getCapon(dataCube,0,0,mainP.dl,sbl,mainP.tavg,[],mainP.dofb,verbose);
    bf_im = abs(bf_im);
    warning('on');
elseif strfind(bf_method, 'IAA-MBSB')
    [IAAImageAmp, IAAImagePow] = getIAAMultiBeam(dataCube,0,0,V, ...
        mainP.P.Rx.Theta,scanGridTheta,pitchInLambdas,1,10,1,verbose,0);
    bf_im = sqrt(IAAImageAmp(:,:,10));
elseif strfind(bf_method, 'IAA-MBMB')
    [IAAImageAmp, IAAImagePow] = getIAAMultiBeam(dataCube,0,0,V, ...
        mainP.P.Rx.Theta,scanGridTheta,pitchInLambdas,1,10,0,verbose,0);
    bf_im = sqrt(IAAImagePow(:,:,10));
else
    error(strcat(bf_method, ' is not a valid beamformer option!'))
end

end

