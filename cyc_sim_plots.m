close all
clc

day = [1, 90, 180, 270, 365];

interface = [-41.04, -39.29, -39.03, -38.90, -38.08;
              -38.76, -36.49, -36.17, -35.98, -35.84;
              -34.78, -31.54, -31.05, -30.77, -30.55];
            
rad_peak = [.0, 0.001, 0.01, 0.012, 0.014;
            10.19, 10.37, 10.37, 10.38, 10.39;
            35.33, 35.26, 35.24, 35.23, 35.24];
          
circ_peak = [207.89, 207.29, 207.20, 207.14, 207.10;
             426.75, 420.23, 419.28, 418.71, 418.29;
             819.74, 800.93, 798.43, 796.97, 795.88];
           
sr_hub = [.271, .256, .255, .253, .253;
          .044, .056, .057, .058, .059;
          .524, .550, .554, .556, .557];
        
sr_rim = [.160, .155, .154, .153, .153;
          .242, .236, .235, .235, .234;
          .760, .740, .737, .735, .734];
          
marker = 'o+v*<s>dx.^p';
line = {'-', '--', '-.'};
for k = 1:3
  figure(1), hold on
  plot(day,(interface(k,:)./interface(k,1)), [marker(k),line{k}], 'Linewidth', 1.5)
  
  
  
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
ylabel('Normalized Interfacial Stress [{\it\sigma_i /\sigma_1}]')
set(gca, 'FontSize', 12)

figure(2), hold on
plot(day, (rad_peak(2,:)./rad_peak(2,1)), [marker(2),line{2}], 'Color', 	[1 0 0], 'Linewidth', 1.5)
plot(day, (rad_peak(3,:)./rad_peak(3,1)), [marker(3),line{3}], 'Color', [0.9290 0.6940 0.1250], 'Linewidth', 1.5)
grid on
legend({'\omega_{Pmin}', '\omega_{Pint}', '\omega_{Pmax}'}, 'Location', 'southeast')
xlabel('Day')
ylabel('Normalized Radial Peak Tensile Str [{\it\sigma_i /\sigma_1}]')
set(gca, 'FontSize', 12)

figure(3), grid on
legend({'\omega_{Pmin}', '\omega_{Pint}', '\omega_{Pmax}'}, 'Location', 'southeast')
xlabel('Day')
ylabel('Normalized Circ. Peak Str [{\it\sigma_i /\sigma_1}]')
set(gca, 'FontSize', 12)

figure(4), grid on
legend({'\omega_{Pmin}', '\omega_{Pint}', '\omega_{Pmax}'}, 'Location', 'southeast')
xlabel('Day')
ylabel('Normalized Hub Peak SR [/]')
set(gca, 'FontSize', 12)

figure(5), grid on
legend({'\omega_{Pmin}', '\omega_{Pint}', '\omega_{Pmax}'}, 'Location', 'southeast')
xlabel('Day')
ylabel('Normalized Rim Peak SR [/]')
set(gca, 'FontSize', 12)







