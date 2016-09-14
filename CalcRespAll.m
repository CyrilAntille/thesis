%
% $Id:$
%
% -----------------------------------------------------------------
%
% Function     : function Data = CalcRasp(P);
%
% Description  : Calculates the respons using Field II
%
% Input        : Struct P, defined in Parameters.m
%
% Language     : Matlab R2008b
%
% Written by   : Andreas Austeng, Ifi-UiO
%
% Version no.  : 1.0  AA 01.07.09 first version
% -----------------------------------------------------------------

function Data = CalcRespAll(P,Phantom)
% if ~isdeployed
%   help CalcRespAll;
%   end

  %% Init Field II. 
  evalc('set_field(''fs'',P.fs)');
  set_field('c',P.c);
  % set_field('use_triangles',1);
  
  % Define the transducers
  Tx = eval(P.Tx.xdc);
  Rx = eval(P.Rx.xdc);
  
  % Impulseresponse and excitation
  xdc_impulse (Tx, P.Tx.ImpulseResp);
  xdc_excitation (Tx, P.Tx.Excitation);
  xdc_impulse (Rx, P.Rx.ImpulseResp);

  % Calculate length of excitation x impulsresponses to be able to do depth
  % correction
  L1 = length(P.Tx.Excitation);
  L2 = length(P.Tx.ImpulseResp);
  L3 = length(P.Rx.ImpulseResp);
  LConvs = -((L1+L2-1) + L3 -1)/2;
  
  % Do calculation
  
  %  Do the calculation
  [Data.image,Data.times]=calc_scat_all (Tx, Rx, Phantom.positions, ...
      Phantom.amplitudes,P.DesFactor);
  
  %   Store the result
  Data.LConvs = LConvs;
  
  
  