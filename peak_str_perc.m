t = [0, 0.5, 1, 5, 10];
inter = [26.98, 19.60, 19.01, 17.66, 17.08];
rad_peak = [50.963, 50.861, 50.831, 50.747, 50.705];
circ_peak = [975.51, 936.54, 933.20, 925.40, 922.01];

al_sr = [.836, .899, .904, .916, .921];
cfrp_sr = [1.01, .975, .972, .965, .961];

dinter = (inter./inter(1));
drad_peak = (rad_peak./rad_peak(1));
dcirc_peak = (circ_peak./circ_peak(1));

dal = (1 - al_sr./al_sr(1))*100;
dcfrp = (1 - cfrp_sr./cfrp_sr(1))*100;

figure(), hold on
plot(t, dinter, 'o-', 'Linewidth', 1.5)
plot(t, drad_peak, '*--', 'Linewidth', 1.5)
plot(t, dcirc_peak, '^-.', 'Linewidth', 1.5)

legend({'Interficial str', 'Peak Radial Str', 'Peak Circ. Str'}, 'Location', 'southwest')
xlabel('Time [y]')
ylabel('Normalized stress [\sigma/\sigma_0]')
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