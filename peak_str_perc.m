close all
t = [0, 0.5, 1, 5, 10];
inter = [32.90, 26.63, 26.13, 24.98, 24.48];
rad_peak = [73.83, 71.36, 71.13, 70.60, 70.36];
circ_peak = [1148.5, 1103.8, 1099.99, 1091.04, 1082.15];

al_sr = [.9253, .9944, 1.000, 1.0128, 1.0182];
cfrp_sr = [.9519, .9189, .9160, .9088, .9057];

dinter = (inter./inter(1));
drad_peak = (rad_peak./rad_peak(1));
dcirc_peak = (circ_peak./circ_peak(1));

dal = (1 - al_sr./al_sr(1))*100;
dcfrp = (1 - cfrp_sr./cfrp_sr(1))*100;

figure(), hold on
plot(t, dinter, 'o-', 'Linewidth', 1.5)
plot(t, drad_peak, '*--', 'Linewidth', 1.5)
plot(t, dcirc_peak, '^-.', 'Linewidth', 1.5)

legend({'Interfacial str', 'Peak Radial Str', 'Peak Circ. Str'}, 'Location', 'southwest')
xlabel('Time [y]')
ylabel('Normalized Peak Stress')
grid on
set(gca, 'Fontsize', 12)

figure(), hold on
plot(t, al_sr, 'r+-', 'Linewidth', 1.5)
plot(t, cfrp_sr, 'kv--', 'Linewidth', 1.5)
legend({'Al Hub', 'CFRP Rim'}, 'Location', 'northwest')
xlabel('Time [y]')
ylabel('Peak SR')
grid on
set(gcs, 'Fontsize', 12)