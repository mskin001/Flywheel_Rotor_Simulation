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
           
sr_hub = [.361, .344, .342, .341, .339;
          .057, .062, .063, .064, .065;
          .689, .719, .724, .726, .723];
        
sr_rim = [.168, .163, .163, .162, .162;
          .252, .25, .249, .249, .249;
          .741, .728, .726, .724, .723];
          
marker = 'o+v*<s>dx.^p';
line = {'-', '--', '-.'};
for k = 1:3
  figure(1), hold on
  plot(day,interface(k,:), [marker(k),line{k}], 'Linewidth', 1.5)
  
  figure(2), hold on
  plot(day, rad_peak(k,:), [marker(k),line{k}], 'Linewidth', 1.5)
  
  figure(3), hold on
  plot(day, circ_peak(k,:), [marker(k),line{k}], 'Linewidth', 1.5)
  
  figure(4), hold on
  plot(day, sr_hub(k,:), [marker(k),line{k}], 'Linewidth', 1.5)
  
  figure(5), hold on
  plot(day, sr_rim(k,:), [marker(k),line{k}], 'Linewidth', 1.5)
end

figure(1), grid on
legend({'\omega_{Pmin}', '\omega_{Pint}', '\omega_{Pmax}'}, 'Location', 'southeast')
xlabel('Day')
ylabel('Interficial Stress [MPa]')
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







