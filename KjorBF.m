function KjorSimMCC(ProsNr);

% Funksjon som startes vha Condor og kj�rer forskjellige BF-metoder
% Variabelen 'ProsNr' er en integer (mellom 1 og "tot. ant. prosesser",
% og denne er identisk til omr�der hvor datafiler ligger


  % M� rette opp p� variabeltype hvis kallet som stand-alone
  if isdeployed
    ProsNr = sscanf(ProsNr,'%d');
  end
% Omraader = [98 82 19 22 39 72 74 76 80 86 89]
% Omraader = [82 76]
 
  R = ProsParams;

  % G�r ned i underomr�de relativ til "Initialdir" satt i Condor "cmd" fil
  % cmd= ['cd ',int2str(Omraader(ProsNr+1))];
  cmd= ['cd ',int2str(ProsNr)];
  eval(cmd)
  
  % Lager "NegAmps" bilder
  DoBeamform('BFDataNegAmps','NegAmps',R);

  % Lager "PosAmps" bilder
  DoBeamform('BFDataPosAmps','PosAmps',R);

