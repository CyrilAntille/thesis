%
% $Id:$
%
% -----------------------------------------------------------------
%
% Function     : function BFImage = BeamformAll(P,Data);
%
% Description  : Do "Delay-and-"-beamforming of data from FieldII's
%                calc_scat_all()
%
% Input        : Struct P, defined in Parameters.m
%                NB: P.Tx.ElPos and P.Rx.ElPos must be defined!!
%                Struct Data defined in CalcRespAll
%
% Language     : Matlab R2008b
%
% Written by   : Andreas Austeng, Ifi-UiO
%
% Version no.  : 1.0  AA 01.07.09 first version
%                1.1  AA 25.10.11 Added P.Tx.Center in calc of xFoc
%                                 Changed '*' to '.*' for xFoc and zFoc
%                                 to be able to have dynamic of multiple foci
%                                 (if P.Tx.FocRad is a vector).
%                                 Added P.Radius as a possibility
%                                 Added P.Tx.Apod as a possibility
%                1.2  AA 22.11.11 Changed P.Tx.Theta to P.Tx.SinTheta
%                                 Added possibility for MLA on Rx
%                                 Added Hilbert of data
% -----------------------------------------------------------------
function BFImage = BeamformAll(P,Data)

% Rearrange data into Radius*Rx*Tx-matrix
RawData = reshape(Data.image,size(Data.image,1), ...
    P.Rx.no_elements,P.Tx.no_elements);

% Take the hilbert-transform of all data! This is easier to get right if
% done on the data before than after BF (if BF-data contains regions with
% zeros).

for ii = 1:P.Tx.no_elements
    RawData(:,:,ii) = hilbert(RawData(:,:,ii));
end


if isfield(P,'Radius'),
    Radius = P.Radius;
    MaxRadius = max(Radius);
    NRadius = length(Radius);
else
    % Define the radiuses to use
    MaxRadius = (Data.times*P.c + Data.LConvs/P.DesFactor/P.fs*P.c + ...
        (size(Data.image,1)-1)/P.fs*P.DesFactor*P.c)/2;
    Radius = [P.MinRadImage:P.c/P.fs*P.DesFactor:MaxRadius].';
    NRadius = length(Radius);
end

% Tabulate "sin" and "cos" of thetas
SinTheta = P.Tx.SinTheta;
CosTheta = cos(asin(P.Tx.SinTheta));


% Allocate memory
if isfield(P.Rx,'DeltaSinTh')
    nMLA = length(P.Rx.DeltaSinTh);
else
    P.Rx.DeltaSinTh = 0;
    nMLA = 1;
end
BFImage.image = zeros(nMLA*length(P.Tx.SinTheta),NRadius,P.Rx.no_elements);
BFImage.Radius = Radius;

% Tx-apodization
if isfield(P.Tx,'Apod'),
    Apod = P.Tx.Apod;
else
    Apod = ones(1,P.Tx.no_elements);
end

% Define CenterPoint for Focus
if isfield(P.Tx,'FocCenter'),
    TFocC = P.Tx.FocCenter;
else
    TFocC = 0;
end
if isfield(P.Rx,'FocCenter'),
    RFocC = P.Rx.FocCenter;
else
    RFocC = 0;
end


% loop over SinThetas
for iTheta = 1:length(P.Tx.SinTheta)
    
    % If P.Tx.Center ~= 0, the tx-radius is wrong:
    if (isfield(P.Tx,'Center') && P.Tx.Center ~= 0),
        % Calculate TxDelays for this SinTheta
        xTx = Radius*SinTheta(iTheta);
        zTx = Radius*CosTheta(iTheta);
        
        % Radius til point that's being imaged
        TxRad = sqrt((xTx - P.Tx.Center).^2 + zTx.^2);
    else
        TxRad = Radius;
    end
    
    % Then loop over all Rx-direction for this Tx-direction
    for iMLA = 1:length(P.Rx.DeltaSinTh)
        % Calculate RxDelay for this direction
        xRx = Radius*(SinTheta(iTheta)+P.Rx.DeltaSinTh(iMLA)) + RFocC;
        zRx = Radius*cos(asin(SinTheta(iTheta)+P.Rx.DeltaSinTh(iMLA)));
        RxDelays = P.fs /P.c / P.DesFactor *...
            sqrt((repmat(xRx,1,P.Rx.no_elements) - ...
            repmat(P.Rx.ElPos(1,:),NRadius,1)).^2 + ...
            (repmat(zRx,1,P.Rx.no_elements) - ...
            repmat(P.Rx.ElPos(3,:),NRadius,1)).^2);
        
        % Loop over Tx
        for iTx = 1:P.Tx.no_elements
            % Find focal-point, fix-focus or dyn-focus
            zFoc = P.Tx.FocRad*CosTheta(iTheta);
            xFoc = P.Tx.FocRad.*SinTheta(iTheta) + TFocC;
            if (isfield(P.Tx,'Center') && P.Tx.Center ~= 0),
                TxFocRad = sqrt((xFoc-P.Tx.Center).^2 + zFoc.^2);
            else
                TxFocRad = P.Tx.FocRad;
            end
            
            % Calculate TxDelays, fix-fokus:
            
            %    TxDelays = repmat(RxDelays(:,iTx),1,P.Tx.no_elements); % Test!
            TxDelays = P.fs /P.c / P.DesFactor *...
                repmat(sqrt((xFoc - P.Tx.ElPos(1,iTx)).^2 + ...
                (zFoc - P.Tx.ElPos(3,iTx)).^2) + ...
                TxRad - TxFocRad,1,P.Rx.no_elements);
            
            TwoWayDelay = RxDelays + TxDelays - ...
                Data.times/P.DesFactor*P.fs-Data.LConvs/P.DesFactor;
            
            for iRx = 1:P.Rx.no_elements
                % keyboard
                Resp = interp1q([1:size(RawData,1)].', ...
                    RawData(:,iRx,iTx),TwoWayDelay(:,iRx)).';
                NotNaN = ~isnan(Resp);
                BFImage.image(nMLA*(iTheta-1)+iMLA,NotNaN,iRx) = ...
                    BFImage.image(nMLA*(iTheta-1)+iMLA,NotNaN,iRx)  ...
                    + Apod(iTx)*Resp(NotNaN);
                
            end
        end
    end
    
%     fprintf('\r Done calculating %.1f%%',100*iTheta/length(P.Tx.SinTheta));
    
end



