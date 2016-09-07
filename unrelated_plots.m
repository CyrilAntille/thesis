%% TEST

mainP = MainParameters();
fprintf('\n0. NTheta: %d.\n', mainP.P.Tx.NTheta)
parfor t=1:4
    fprintf('\n1. %d: NTheta: %d.\n', t, mainP.P.Tx.NTheta)
    mainP.set.P(101*t);
    fprintf('\n2. %d: NTheta: %d.\n', t, mainP.P.Tx.NTheta)
end

%% Aperture smoothing function

M = 5; % Number of sensors
d = 1; % Distance between each sensor, in meters
Dm = 1; % Size of sensor array, in meters

kx = linspace(-3, 3, 100); % Wavenumber, in pi/d
k = kx.*(pi/d);

Wa = sin(k.*(M*d/2))./sin(k.*(d/2));
We = sin(k.*(Dm/2))./(k./2);
Wt = We .* Wa;

figure;
plot(kx, Wa, 'b', 'LineWidth', 2, 'LineStyle', '-.')
hold on
plot(kx, We, 'r', 'LineWidth', 2, 'LineStyle', '--')
plot(kx, Wt, 'k', 'LineWidth', 3, 'LineStyle', '-')

legend({'W(k_x) array', 'W(k_x) sensor', 'W(k_x) total'}, 'Location', 'best')
xlabel('Wavenumber k_x in [\pi/d]')
ylabel('Amplitude')

hold off;

figure;
max_db = max(db(Wt));
plot(kx, db(Wa)-max(db(Wa)), 'b', 'LineWidth', 2, 'LineStyle', '-.')
hold on
plot(kx, db(We)-max(db(We)), 'r', 'LineWidth', 2, 'LineStyle', '--')
plot(kx, db(Wt)-max(db(Wt)), 'k', 'LineWidth', 3, 'LineStyle', '-')

legend({'W(k_x) array', 'W(k_x) sensor', 'W(k_x) total'}, 'Location', 'best')
xlabel('Wavenumber k_x in [\pi/d]')
ylabel('Power [dB]')
ylim([-55 0])

hold off; pause; close all

%%
figure; hold on;
for m=1:length(methods_set)
    m_BF = data_BF{m};
    for b=1:length(num_beams)
        Pb = copyStruct(P);
        Pb.Tx.NTheta = num_beams(b);
        Pb.Tx.SinTheta = linspace(-Pb.Tx.SinThMax,Pb.Tx.SinThMax,Pb.Tx.NTheta);
        Pb.Tx.Theta = asin(Pb.Tx.SinTheta);

        b_BF = m_BF{b};
        b_DA = data_DA{b};
        
        bf_img = db(b_BF)  - max_peak_DAS;
        z_sep = find(b_DA.Radius >= P.Tx.FocRad, 1);
        beampattern_m = bf_img(z_sep,:);
        plot(Pb.Tx.Theta, beampattern_m)
    end
end
hold off; legend(methods_set); pause; close;

