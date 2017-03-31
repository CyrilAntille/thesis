%
% $Id:$
%
% -----------------------------------------------------------------
%
% Function     : function P = Parameters;
%
% Description  : Defines the parameters used for data-generation
%
% Language     : Matlab R2008b
%
% Written by   : Andreas Austeng, Ifi-UiO
%
% Version no.  : 1.0  AA 01.07.09 first version
% -----------------------------------------------------------------

function P = Parameters(NTheta)
%   help Parameters;
    if nargin == 0
        NTheta = 161;
    end

  %%
  % Medium parameters
  
  P.c     =  1500;    % Speed-of-sound [m/s]
  P.fc    = 3e6;     % Transducer center frequency [Hz]
  lambda = P.c/P.fc; % -> 0.5 mm
  %%
  % Field II parametres
  P.fs    =  90e6;   % Sampling frequency user by FieldII [Hz]
  P.DesFactor = 4;    % Decimation factor
%   P.DesFactor = 1;    % Decimation factor
  
  %%
  % Parameters defining Tx and Rx

  % Tx:
   P.Tx.xdc = ['xdc_linear_array(P.Tx.no_elements, P.Tx.width, ' ...
      'P.Tx.height, P.Tx.kerf, P.Tx.no_sub_x, ' ...
      'P.Tx.no_sub_y, P.Tx.focus);'];   % The function defining the aperture


  % Must provide all input-data for this function
  P.Tx.no_elements = 96;                    % No of physical elements in x
  P.Tx.height      = 10e-3;                 % Height of elements
  P.Tx.pitch       = lambda/2;                % Center-to-center distance 
  P.Tx.kerf        = 1e-5;                  % Distance between elements  
  P.Tx.width       = P.Tx.pitch-P.Tx.kerf;  % Width of element in x
  P.Tx.no_sub_x    = 3;                     % No of sub-division in x 
  P.Tx.no_sub_y    = 8;                    % No of sub-division in x 
  P.Tx.focus       = [0 0 40]*1e-3;         % Initial electonic focus
  P.Tx.Foc = 2;

  P.Tx.ElPos = [[-(P.Tx.no_elements-1)/2:(P.Tx.no_elements-1)/2]*P.Tx.pitch; ...
		zeros(2,P.Tx.no_elements)];

  % Rx:
  P.Rx.xdc = ['xdc_linear_array(P.Rx.no_elements, P.Rx.width, ' ...
      'P.Rx.height, P.Rx.kerf, P.Rx.no_sub_x, ' ...
      'P.Rx.no_sub_y, P.Rx.focus);'];   % The function defining the aperture
  

      % Must provide all input-data for this function
  P.Rx.no_elements = 96;                    % No of physical elements in x
  P.Rx.height      = 10e-3;                 % Height of elements
%   P.Rx.pitch       = P.c/P.fc /2;           % Center-to-center distance 
%   P.Rx.pitch       = 0.33e-3;                % Center-to-center distance 
  P.Rx.pitch       = lambda/2;                % Center-to-center distance 
  P.Rx.kerf        = 1e-5;                  % Distance between elements  
  P.Rx.width       = P.Rx.pitch-P.Rx.kerf;  % Width of element in x
  P.Rx.no_sub_x    = 3;                     % No of sub-division in x 
  P.Rx.no_sub_y    = 8;                    % No of sub-division in x 
  P.Rx.focus       = [0 0 40]*1e-3;         % Initial electonic focus
  
  P.Rx.ElPos = [[-(P.Rx.no_elements-1)/2:(P.Rx.no_elements-1)/2]*P.Rx.pitch; ...
		zeros(2,P.Rx.no_elements)];

% Impulsresponse, Tx and Rx
% 3-periods and Blackman gives ~78% BW
% 3-periods and Hanning  gives ~65.5% BW
% 3-periods and Hamming  gives ~61 % BW
% 2-periods and Hamming  gives ~93% BW

  ImpulseResp = -sin(2*pi*P.fc*(0:1/P.fs:(3-1/P.fc)/P.fc));
  P.Tx.ImpulseResp = ImpulseResp.*blackman(max(size(ImpulseResp))).';
  % P.Tx.ImpulseResp = 1;
  
  P.Rx.ImpulseResp = ImpulseResp.*blackman(max(size(ImpulseResp))).';
  % P.Rx.ImpulseResp = 1;
  
  % Excitation
  %   P.Tx.Excitation = sin(2*pi*P.fc*(0:1/P.fs:3/P.fc));
  P.Tx.Excitation = square(2*pi*P.fc*(0:1/P.fs:(1-1/P.fs)/P.fc));
  % P.Tx.Excitation = 1;
  
  %%
  % Scaling of impulseresponse and excitation to be able to use single
  % precision. OBS: This must be checked for different array setups!
  
  P.Tx.Excitation = P.Tx.Excitation * 1e6;
  P.Tx.ImpulseResp = P.Tx.ImpulseResp/max(P.Tx.ImpulseResp) * 1e6;
  P.Rx.ImpulseResp = P.Rx.ImpulseResp/max(P.Rx.ImpulseResp) * 1e6;
  
  
  %% 
  % Define imaging sector and focus. Theta uniform in sin(Theta)!
  P.Tx.NTheta = NTheta;
  P.Tx.SinThMax = 0.3; 
  P.Tx.SinTheta = linspace(-P.Tx.SinThMax,P.Tx.SinThMax,P.Tx.NTheta);
  P.Tx.Theta = asin(P.Tx.SinTheta);
  
  
  P.Tx.FocRad   = 40e-3;    % Focal radius
  P.Rx.NoMLA    = 1;        % No of Rx-lines per Tx-line 
  P.Rx.DeltaTh  = 0;        % Rx-offset relative to Tx-line [rad]
  P.MinRadImage = 35e-3;     % Minimum radius for which data is recorded [m]


