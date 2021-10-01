rim = [.16 .2 .33];
w_max = 24150 * pi / 30;
w_min = 6037.5 * pi / 30;
p = [0, 100, 250, 500, 750, 1000, 1250, 1500, 1750, 2000];%, 5000, 10000,...
    %25000, 75000, 100000, 200000, 300000, 400000, 500000];% 600000,...
    %700000, 800000, 900000, 1000000];
sh1 = -1 * [0, 5.72E-04, 1.57E-03, 2.86E-03, 4.29E-03, 5.71E-03, 7.14E-03,...
     8.58E-03, 1.00E-02, 1.14E-02];%, 2.86E-02, 5.72E-02, 1.43E-01, 2.87E-01,...
     %4.30E-01, 5.75E-01, 1.16E+00, 1.75E+00, 2.34E+00, 2.95E+00];%, 3.56E+00,...
     %4.18E+00, 4.81E+00, 5.45E+00, 6.10E+00]; in [MPa]
     
sh2 = [0, 2.28e-3, 5.72e-3, 1.14e-2, 1.72e-2, 2.29e-2, 2.86e-2, 3.43e-2,...
    4.00e-2, 4.57e-2];
rho = [2810, 1560];
h = .43;

for k = 1:length(rim) - 1
    in(k) = 0.5 * rho(k) * pi * h * (rim(k+1)^4 - rim(k)^4);
end
I = sum(in);

dec = -p / (I * w);
acc = p / (I * w);

yyaxis left
hold on
shr_dec = plot(-p, sh1*10^3, 'b-o', 'Linewidth', 1.5);
shr_acc = plot(p, sh2*10^3, 'b--*', 'Linewidth', 1.5);
ylabel('Peak Shear Str [kPa]')

yyaxis right
hold on
pwr_dec = plot(-p, dec*(60^2/(2*pi)), 'r-v', 'Linewidth', 1.5);
pwr_acce = plot(p, acc*(60^2/(2*pi)), 'r-v', 'Linewidth', 1.5);
legend([shr_dec, shr_acc, pwr_dec] , '24,150 rpm', '6037.5 rpm',...
    'acceleration')
ylabel('Acceleration [rad/sec^2]')
xlabel('Absolute Power [kW]')
ylim([-200 200])

grid on
set(gca, 'XLim', [-2000,2000], 'FontSize', 12)