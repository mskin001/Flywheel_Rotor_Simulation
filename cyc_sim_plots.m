close all
clc

day = [1, 90, 180, 270, 365];

interface = [-40.04, -38.47, -38.27, -38.16, -38.07;
              -37.26, -35.23, -34.95, -34.78, -34.66;
              -32.36, -29.44, -29.0, -28.74, -28.56];
            
rad_peak = [.797, .914, .931, .941, .948;
            19.49, 19.32, 19.29, 19.27, 19.26;
            57.58, 56.58, 56.41, 56.31, 56.24];
          
circ_peak = [211.61, 211.7, 211.69, 211.69, 211.68;
             475.45, 467.96, 466.87, 466.22, 465.74;
             945.41, 926.11, 923.07, 921.29, 919.96];
           
sr_hub = [.361, .3441, .342, .3406, .3396;
          .057, .0622, .0633, .064, .0645;
          .686, .7187, .7236, .7264, .7285];
        
sr_rim = [.168, .1634, .1628, .1624, .1622;
          .2519, .2498, .2494, .2492, .2491;
          .7411, .7278, .7257, .7243, .7234];
          
marker = 'o+v*<s>dx.^p';
line = {'-', '--', '-.'};
for k = 1:3
  figure(1), hold on
  plot(day,(interface(k,:)./interface(k,1)), [marker(k),line{k}], 'Linewidth', 1.5)
  
  figure(2), hold on
  plot(day, (rad_peak(k,:)./rad_peak(k,1)), [marker(k),line{k}], 'Linewidth', 1.5)
  
  figure(3), hold on
  plot(day, (circ_peak(k,:)./circ_peak(k,1)), [marker(k),line{k}], 'Linewidth', 1.5)
  
  figure(4), hold on
  plot(day, (sr_hub(k,:)./sr_hub(k,1)), [marker(k),line{k}], 'Linewidth', 1.5)
  
  figure(5), hold on
  plot(day, (sr_rim(k,:)./sr_rim(k,1)), [marker(k),line{k}], 'Linewidth', 1.5)
end

figure(1), grid on
legend({'\omega_{Pmin}', '\omega_{Pint}', '\omega_{Pmax}'}, 'Location', 'southeast')
xlabel('Day')
ylabel('Interfacial Stress')
set(gca, 'FontSize', 12)

figure(2), grid on
legend({'\omega_{Pmin}', '\omega_{Pint}', '\omega_{Pmax}'}, 'Location', 'southeast')
xlabel('Day')
ylabel('Radial Peak Stress [MPa]')
set(gca, 'FontSize', 12)

figure(3), grid on
legend({'\omega_{Pmin}', '\omega_{Pint}', '\omega_{Pmax}'}, 'Location', 'southeast')
xlabel('Day')
ylabel('Circ. Peak Stress [MPa]')
set(gca, 'FontSize', 12)

figure(4), grid on
legend({'\omega_{Pmin}', '\omega_{Pint}', '\omega_{Pmax}'}, 'Location', 'southeast')
xlabel('Day')
ylabel('Hub Peak SR')
set(gca, 'FontSize', 12)

figure(5), grid on
legend({'\omega_{Pmin}', '\omega_{Pint}', '\omega_{Pmax}'}, 'Location', 'southeast')
xlabel('Day')
ylabel('Rim Peak SR')
set(gca, 'FontSize', 12)







