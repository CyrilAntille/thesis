function [Data] = LagCystPhantom(ProsNr,fname,VesselAmpl,NPts,Seed, ...
    SaveToFile, P, mmShift);

% Funksjon som redefinerer Tx-aperturen til ? v?re et element, 
% kaller deretter PlanewaveVesselPhantom og CalcScatAll for dette ene 
% elementet og lagrer utfiler
%
% Funkjonen er tenkt kalt fra Condor-srcipt, og antall condor-prosesser m?
% v?re lik antall elementer i Tx-apertruren.

%addpath '/hom/dsb/projects/matlab/beamforming/functions/BeamForm/' -end
  
   % Tolker input-variable rett
    if isdeployed,
      ProsNr = sscanf(ProsNr,'%d');
      VesselAmpl = sscanf(VesselAmpl,'%f'); 
      NPts = sscanf(NPts,'%f');
      Seed = sscanf(Seed,'%f');
    end

    
    %field_init(0);


%P = Parameters;
% ElPos = P.ElPos;
% NumEls = P.Tx.no_elements;

% F?rst genererer Tx-aperture vha field, deretter lese ut data for et og et
% element og lage ny beskrivese av Tx-aperture inneholdende denne infoen.

%% Init Field II. 
% set_field('fs',P.fs);
% set_field('c',P.c);
% % set_field('use_triangles',1);
%   
% % Define the transducers og leser ut info om denne
% Tx_tmp = eval(P.Tx.xdc);
% Txdata=xdc_get(Tx_tmp,'rect').';
% xdc_free(Tx_tmp);
% 
% % Lager ny definisjon av 'xdc'
% P.Tx.xdc_orig = P.Tx.xdc;  % orginal xdc
% P.Tx.xdc = ['xdc_rectangles(P.Tx.rect,P.Tx.centers,P.Tx.focus);']
% 
% 
% % Definerer 'rect' og 'center' til bruk i xdc_rectangles for et element
% RowsToGet = [11:22 5 3:4 8:10];  
% 
% rect = [];
% centers = [];
% TheseEls = [];
% ENR = 0;
% 
% for ii = 2*ProsNr+[0 1]
%     
%     ENR = ENR + 1;
%     
%     ElNr = find(Txdata(:,1)==ii);
%     
%     
%     rect =  [rect; ENR*ones(length(ElNr),1) Txdata(ElNr,RowsToGet)];
%     centers = [centers; Txdata(ElNr(1),24:26)];
%     
%     
%     TheseEls = [TheseEls ii]; % Lagrer info om hvilket element dette egentlig er.
% 
% end
% 
% P.Tx.rect = rect;
% P.Tx.centers = centers;
% P.Tx.TheseEls = TheseEls;

% if 0,
[Phantom] = PlanewaveVesselPhantom(P,VesselAmpl,NPts,Seed);
% else
%   points = [0 0 40]*1e-3;
%   Phantom = PointPhantom(points);
% end

field_init(0);
Data = CalcRespAll(P,Phantom);

% field_end();
if SaveToFile
    tmp = ['save -v7.3 ',fname,'_Tx',int2str(ProsNr+1),...
           ' ProsNr fname  VesselAmpl NPts Seed P Data;  '];

    eval(tmp);
end



