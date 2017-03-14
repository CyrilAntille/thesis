%
% Generate 2 signals in noise according to page 73 of 
% H. Krim, M. Viberg, "Two decades of array signal processing research
% - The parametric approach," IEEE Signal Processing Magazine, pp.67-94,
% July 1996.
%
% Note that definition of direction of arrival (DOA) is different from that
% of Krim and Viberg in that 0 degrees is normal incidence
%
% Written by:
% S. Holm, Department of Informatics, University of Oslo
% 31. Oct 1997 First version 
% 21. Nov 1997 Corrected for noise in imaginary part 21. November 1997
% 12. Dec 1997 Added random phase noise to signals to better ensure incoherence


M = 100; 				% number of sensors
N = 10* M; 				% number of samples

SNR1 = 0; 				% dB
SNR2 = 0; 				% dB
theta1 = 1; 				% degrees
theta2 = -1; 				% degrees
gamma  = 45;				% degrees, phase shift between sources
d = 0.5; 				% element spacing in array
k = 2*pi; 				% normalized wavenumber
T = 0.5; 				% sampling interval
omega1 = 2*pi; 				% normalized frequency
omega2 = 2*pi;
%*(1-3/N); 			        % one way to make incoherent 
                                        % is to make the frequencies 
                                        % slightly different


randn('seed',0); 			% make sure the same noise
                                        % sequence is generated each time
noise = randn(M,N) + j*randn(M,N);      % Gaussian noise, variance = 1.0
                                        % for real and imag part, i.e.
					% variance of noise = 2
amp1 = sqrt(2)*10^(SNR1/20);       
amp2 = sqrt(2)*10^(SNR2/20);

phi1 = -k*d*sin(theta1*pi/180);
phi2 = -k*d*sin(theta2*pi/180);
a1 = exp(j*phi1).^[0:M-1]';
a2 = exp(j*phi2).^[0:M-1]';
A = [a1 a2];

phase1 = j*2*pi*rand(N,1); % random phase noise to make sources incoherent
phase2 = j*2*pi*rand(N,1); % random phase noise to make sources incoherent
%
% uncomment this line to make sources coherent:
phase2 = phase1;

signal1 = amp1*exp(phase1).*exp(j*omega1*T).^[0:N-1]';
signal2 = amp2*exp(j*gamma*pi/180)*exp(phase2).*exp(j*omega2*T).^[0:N-1]';
s = [signal1 signal2]';

x = A*s + noise;

disp(['Signal 1:       amplitude ',num2str(amp1),', frequency ', ...
	num2str(omega1/(2*pi)),', direction of arrival ', ...
	num2str(theta1),' deg'])
disp(['Signal 2:       amplitude ',num2str(amp2),', frequency ', ...
	num2str(omega2/(2*pi)),', direction of arrival ', ...
	num2str(theta2),' deg, rel. phase ',num2str(gamma),' deg'])
disp(['Noise: standard deviation ',num2str(sqrt(1+1))])
disp(['Output vector x is size ',num2str(M),' elements by ', ...
    num2str(N),' time samples'])


