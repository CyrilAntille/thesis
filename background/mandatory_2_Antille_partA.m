% INF5410 - Mandatory exercise 2
% Cyril Antille
%% Part A: High Resolution Beamforming
%Generate input data (variable x)
proj2b

%%
L_set = [M, 3*M/4, M/2];
R_set = cell(size(L_set));
for lidx = 1:length(L_set)
    L = L_set(lidx); % subarrays length
    numK = M - L + 1; % number of subarrays
    Rx = zeros(L);
    %Rx = x*x';
    for n=1:N % Temporal averaging
        for kidx = 1:numK % spatial averaging
            Rx = Rx + x(kidx:kidx+L-1,n) * x(kidx:kidx+L-1,n)';
        end
    end
    Rx = Rx ./ (numK * N);
    R_set{lidx} = Rx;
%     imagesc(real(Rx)); colorbar
%     xlabel('Sensor Number'); ylabel('Sensor Number');
%     title(strcat('Spatial Covariance Matrix. L =  ', int2str(L)))
%     pause;
end
% close;

%% 1) temporal estimate of spatial correlation matrix (slide 77)

% 2,3,5) Spatial spectrum
theta_range = linspace(-5,5,200); %degrees

% Building steering vector a
% Following Prof. Holm's definition of DOA -> change cos function from page
% 73 to sin function: cos(theta + 90) = - sin(theta)
phi_range = k * d * sind(theta_range);
a_matrix = zeros(M,length(theta_range));
for m_index = 1:M
    a_matrix(m_index,:) = exp(1i*(m_index-1)*phi_range);
end

% DAS beamformer
classical_spectrum = zeros(1,length(theta_range));
Rx = R_set{1}; % Assumes first Rx is within spatial averaging
for theta_index = 1:length(theta_range)
    a = a_matrix(:,theta_index);
    % 2) Conventional method
    classical_spectrum(theta_index) = (a' * Rx * a) ./ (a' * a);
end
norm_factor = max(classical_spectrum);

lgd = {strcat('DAS method. M=', int2str(M))};
linestyle_list = {'-.','--','-',':', '-'};
markers_list = {'+','x','diamond','o'};
color_list = {'b','r','g','k','m','c'};
figure;
plot(theta_range, db(real(classical_spectrum./norm_factor)), 'Color', ...
    color_list{1}, 'LineWidth', 1.5, 'LineStyle', linestyle_list{1});
% 'Marker', markers_list{1});
hold on;

for lidx = 1:length(L_set)
    L = L_set(lidx); % subarrays length
    Rx = R_set{lidx};
    spectrum = zeros(1,length(theta_range));
    for theta_index = 1:length(theta_range)
        a = a_matrix(1:L,theta_index);
        spectrum(theta_index) = (a' * a) / (a' * inv(Rx) * a);
    end
    plot(theta_range, db(real(spectrum./norm_factor)), 'Color', ...
        color_list{lidx+1}, 'LineWidth', 1.5, 'LineStyle', linestyle_list{lidx+1});
%     'Marker', markers_list{lidx+1});
    lgd{end+1} = strcat('Capon L=', int2str(L));
end
 
% TODO: Increase M (and keep N = 100*M), correct Rx, Normalize by DAS max
% For choice of variables in Methods, find citation in INF5410 lecture for
% requirements for robust Rx estimation (time samples >= 100*num sensors).
% Also check Capon 1970

title('Spatial Spectrum Estimates');
xlabel('DOA [deg]'); ylabel('Spectrum [dB]');
xlim([-5 5])
legend(lgd, 'Location', 'best')
pause; close
