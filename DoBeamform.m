function DoBeamform(DataFil,PartOutfileName,R);
  
  % Laster inn data
  cmd = ['load ',DataFil];
  eval(cmd);
  
  % Korrigerer vinkler og dybder
  if ~isempty(R.ProcRadii),
      P.Radii = BFData.Radius(R.ProcRadii);
  else
      P.Radii = BFData.Radius;
      R.ProcRadii = 1:length(BFData.Radius);
  end
  if ~isempty(R.ProcVinkler),
      P.Tx.Theta = P.Tx.Theta(R.ProcVinkler);
  else
      R.ProcVinkler = 1:length(P.Tx.Theta);
  end
  
  % Hilbert-transformerer data og snu på rekkefølge av dimensjoner
  
  [m,n,o]  = size(BFData.image);
  dataCube = complex(zeros(n,m,o));
  for index = 1:size(BFData.image,3)
    dataCube(:,:,index) = hilbert(transpose(BFData.image(:,:,index)));
  end
  
  
  
  
  
  %% Capon-BF
  
  if R.BF_CAPON==1,
    
    for iL = R.MV.L
      % for alle subarrayer
      
      for idelta = R.MV.delta
	% for alle reg. koeffisienter
	
	for iK = R.MV.K
	  % for alle tidsmidlinger
	  
	  %% Kaller capon-funksjon skrevet av Are 
	  [imAmplitude imPower] = getCapon(dataCube, R.ProcRadii, R.ProcVinkler, ...
						     idelta, iL, iK, [], R.FB,1);
	  
	  % Korrigere størrelse på bilder
	  imAmplitude = imAmplitude(R.ProcRadii,R.ProcVinkler);
	  imPower = imPower(R.ProcRadii,R.ProcVinkler);
	  
	  % Lagrer bildet
	  SaveImage('MV', PartOutfileName ,imAmplitude, imPower, P,iL, iK, idelta);
	  
	end
      end
    end
  end
  
  

  
  if R.BF_DAS == 1,
    
    % Gjøre standard sum over mottakselementer
    imAmplitude = abs(sum(dataCube,3))/P.Rx.no_elements;
    
    % Korrigere størrelse på bilder
    imAmplitude = imAmplitude(R.ProcRadii,R.ProcVinkler);
    imPower = [];
    
    % Lagrer bildet
    SaveImage('DAS', PartOutfileName ,imAmplitude, imPower, P, P.Rx.no_elements, 0, 0);
  end
  
  if R.BF_DAS_HAM == 1,
    
    % Gjøre Vektet-sum over mottakselementer
    imAmplitude = zeros(size(dataCube,1), size(dataCube,2));
    Vindu = hamming(P.Rx.no_elements);
    Vindu = Vindu/sum(Vindu);
    for index = 1:size(dataCube,1)
      imAmplitude(index,:) = abs(squeeze(dataCube(index,:,:))*Vindu); 
    end
    
    % Korrigere størrelse på bilder
    imAmplitude = imAmplitude(R.ProcRadii,R.ProcVinkler);
    imPower = [];
    
    % Lagrer bildet
    SaveImage('DASHAM', PartOutfileName ,imAmplitude, imPower, P, P.Rx.no_elements, 0, 0);
  end

  
   %% Apes-BF
  
  if R.BF_APES==1,
    
    for iL = R.Apes.L
      % for alle subarrayer
      
      for idelta = R.Apes.delta
	% for alle reg. koeffisienter
	
	for iK = R.Apes.K
	  % for alle tidsmidlinger
	  
	  %% Kaller capon-funksjon skrevet av Are 
	  
	  [imAmplitude] = getApes(dataCube, R.ProcRadii, R.ProcVinkler, ...
						     idelta, iL, iK, 1, [], R.FB,1);
	  
	  % Korrigere størrelse på bilder
	  imAmplitude = imAmplitude(R.ProcRadii,R.ProcVinkler);
	  imPower = [];
	  % Lagrer bildet
	  SaveImage('Apes', PartOutfileName ,imAmplitude, imPower, P,iL, iK, idelta);
	  
	end
      end
    end
  end
 
  