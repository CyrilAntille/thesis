function R = ProsParams();
% Funksjon som setter opp parametre for prosessering av 
% Theta-Radii-Rx-datakuber for å generere bilder 
% vha forskjellige stråleformede


% Spesifisering av indexer for hvilke radiuser 
R.ProcRadii = 1:332;

% og vinkler som skal prosesseres
 R.ProcVinkler = 1:241;

% Spesifisering av hvilke bilder og hvilke parametre som skal kjøres
R.BF_DAS = 1;
R.BF_DAS_HAM = 1;
% R.BF_CAPON = 1;
R.BF_APES = 1;

 R.BF_DAS = 0;
 R.BF_DAS_HAM = 0;
 R.BF_CAPON = 0;
% R.BF.APES = 0;

R.MV.L = [60 64];
R.MV.delta = [1/100 1/10];
R.MV.K = [0 20];

R.Apes.L = [64];
R.Apes.delta = [1/100];
R.Apes.K = [0];
