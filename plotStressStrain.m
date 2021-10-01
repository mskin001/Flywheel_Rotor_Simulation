function plotStressStrain(legTxt, unit)
%% -----------------------------------------------------------------------------
% Define global variables, arrays, and structures
% ------------------------------------------------------------------------------
global rim rArr plotWhat results mat
uArr = results.uArr;
sArr = results.sArr;
tau =  results.tauArr;

marker = 'o+v*<s>dp^x'; %o+v*<

if strcmp(unit, 'ips')
  lng = 39.3701; %convert length to inches
  lng_unit = 'in';
  force = 0.000145038; %convert pa to psi
  force_unit = 'psi';
elseif strcmp(unit, 'mmMPas')
  lng = 1000;
  lng_unit = 'mm';
  force = 10^-6;
  force_unit = 'MPa';  
else
  lng = 1;
  lng_unit = 'm';
  force = 10^-6;
  force_unit = 'MPa';
end

if strcmp(legTxt{1},'auto')
    time = results.time(1:plotWhat.interval:end);
    time(end+1) = results.time(end);
    for k = 1:length(time)
        time_str{k} = num2str(time(k));
        legTxt{k} = [time_str{k}, ' sec'];
    end
end

%% -----------------------------------------------------------------------------
% Define rim origional centers and radii
% ------------------------------------------------------------------------------
numRims = length(rim);
oc = zeros(numRims,2);
or = rim;

%% -----------------------------------------------------------------------------
% Create a custom plot
% ------------------------------------------------------------------------------
% To create additional custom plots copy and paste this section to create
% as many custom plots as desired.
halfWay = round(length(sArr)/2);

if strcmp(plotWhat.custom1, 'yes')
  figure(), hold on
  plot(rArr*lng, sArr{1}(3,:,1)*force*10,'-', 'Linewidth', 1.5)
  rad_data = csvread('Apa2011_10xradial.csv');
  plot(rad_data(:,1), rad_data(:,2), 'kd', 'MarkerFaceColor', 'k')
  
%   plot(rArr*lng, sArr{2}(3,:,1)*force,'--', 'Linewidth', 1.5)
%   rad_data = csvread('radial-stress_t2_pressfit_cylinder.csv');
%   plot(rad_data(:,1), rad_data(:,2), 'ks', 'MarkerFaceColor', 'k')
%   
%   plot(rArr*lng, sArr{3}(3,:,1)*force,':', 'Linewidth', 1.5)
%   rad_data = csvread('radial-stress_t3_pressfit_cylinder.csv');
%   plot(rad_data(:,1), rad_data(:,2), 'k^', 'MarkerFaceColor', 'k')
  
%   grid on
%   xlabel(['Radius [', lng_unit, ']'])
%   ylabel(['Stress [', force_unit,']'])
%   legend('Model t1', 'Tzeng t1', 'Model t10', 'Tzeng t10', 'Model tinf', 'Tzeng tinf'...
%       ,'Location', 'southeast')
%   set(gca, 'FontSize', 12)
%   
%   figure(), hold on
  plot(rArr*lng, sArr{1}(1,:,1)*force, '-', 'Linewidth', 1.5)
  hoop_data = csvread('Apa2011_hoop.csv');
  plot(hoop_data(:,1), hoop_data(:,2), 'kd', 'MarkerFaceColor', 'k')
    
%   plot(rArr*lng, sArr{2}(1,:,1)*force, '--', 'Linewidth', 1.5)
%   hoop_data = csvread('hoop-stress_t2_pressfit_cylinder.csv');
%   plot(hoop_data(:,1), hoop_data(:,2), 'ks', 'MarkerFaceColor', 'k')
%   
%   plot(rArr*lng, sArr{3}(1,:,1)*force, ':', 'Linewidth', 1.5)
%   hoop_data = csvread('hoop-stress_t3_pressfit_cylinder.csv');
%   plot(hoop_data(:,1), hoop_data(:,2), 'k^', 'MarkerFaceColor', 'k')
  
  plot(rArr*lng, tau{1}*force, ':', 'Color', [0.4940 0.1840 0.5560], 'Linewidth', 1.5)
  tau_data = csvread('aparicio2011_results.csv');
  plot(tau_data(:,1), tau_data(:,2), 'ko')
  grid on
  legend('Model t1', 'Tzeng t1', 'Model t10', 'Tzeng t10', 'Model tinf', 'Tzeng tinf'...
      ,'Location', 'southeast')
  xlabel(['Radius [', lng_unit, ']'])
  ylabel(['Stress [', force_unit,']'])
%   legend('Model 10*\sigma_r', 'Aparicio2011, 10*\sigma_r','Model \sigma_\theta',...
%       'Aparicio2011, \sigma_\theta', 'Model \tau_r_\theta', 'Aparicio2011, \tau_r_\theta',...
%       'NumColumns', 3, 'Location', 'southoutside')
  set(gca, 'FontSize', 12)
  fprintf('Custom plot 1: Complete\n')
end

%% -----------------------------------------------------------------------------
% Make figures
% ------------------------------------------------------------------------------
% -------------- Radial displacements ------------------------------------------
if strcmp(plotWhat.radDis, 'yes')
  try
    subSet = uArr(1:plotWhat.interval:end);
  catch
    subSet = uArr;
  end
  
  radDis = figure('Visible','on');
  hold on
%   plot(rArr*convert, uArr{1}(1,:)*convert, 'LineWidth', 1.5)
  for k = 1:length(subSet)
    plot(rArr*lng, subSet{k}(1,:)*lng, 'LineWidth', 1.5)
  end
  plot(rArr*convert, uArr{end}(1,:)*convert, 'LineWidth', 1.5)
  xlabel(['Radius [', lng_unit, ']'])
  ylabel(['Radial Displacement [', lng_unit, ']'])
  legend(legTxt, 'Location', 'southeast')
  set(gca, 'FontSize', 12)
  grid on
  fprintf('Radial Displacement Plot: Complete\n')
end

% -------------- Radial stress -------------------------------------------------
if strcmp(plotWhat.radStr, 'yes')
  radStr = figure('Visible','on');
  try
    subSet = sArr(1:plotWhat.interval:end);
  catch
    subSet = sArr;
  end

  hold on
  for k = 1:length(subSet)
    plot(rArr*lng, subSet{k}(3,:,1)*force, [marker(k),'-'], 'MarkerIndices', 1:5:length(rArr), 'LineWidth', 1.5);
  end
  plot(rArr*lng, sArr{end}(3,:,1)*force, [marker(k+1),'-'], 'MarkerIndices', 1:5:length(rArr), 'LineWidth', 1.5);
  xlabel(['Radius [', lng_unit, ']'])
  ylabel(['Radial Stress [', force_unit, ']'])
  legend(legTxt, 'Location', 'southeast')
  set(gca,'Ytick', -50:25:50, 'YtickLabel', -50:25:50, 'FontSize', 12)
  grid on
  fprintf('Radial Sress Plot: Complete\n')
end

% -------------- Hoop stress ---------------------------------------------------
if strcmp(plotWhat.hoopStr, 'yes')
  hoopStr = figure('Visible','on'); %#ok<*NASGU>
  try
    subSet = sArr(1:plotWhat.interval:end);
  catch
    subSet = sArr;
  end

  hold on
  for k = 1:length(subSet)
    plot(rArr*lng, subSet{k}(1,:,1)*force, [marker(k),'-'], 'MarkerIndices', 1:5:length(rArr),  'LineWidth', 1.5);
  end
  plot(rArr*lng, sArr{end}(1,:,1)*force, [marker(k+1),'-'], 'MarkerIndices', 1:5:length(rArr),  'LineWidth', 1.5);

  xlabel(['Radius [', lng_unit, ']'])
  ylabel(['Circumferential Stress [', force_unit, ']'])
  legend(legTxt, 'Location', 'southeast', 'NumColumns', 2)
  set(gca, 'FontSize', 12)
  grid on
  fprintf('Hoop Stress Plot: Complete\n')
end

% -------------- Axial stress -------------------------------------------------
if strcmp(plotWhat.axialStr, 'yes')
  radStr = figure('Visible','on');
  try
    subSet = sArr(1:plotWhat.interval:end);
  catch
    subSet = sArr;
  end

  hold on
  for k = 1:length(subSet)
    plot(rArr*lng, subSet{k}(2,:,1)*force, [marker(k),'-'], 'MarkerIndices', 1:5:length(rArr), 'LineWidth', 1.5);
  end
  plot(rArr*lng, sArr{end}(2,:,1)*force, [marker(k+1),'-'], 'MarkerIndices', 1:5:length(rArr), 'LineWidth', 1.5);
  
  xlabel(['Axial [', lng_unit, ']']);
  ylabel(['Axial Stress [', force_unit, ']']);
  legend(legTxt, 'Location', 'southeast')
  set(gca, 'FontSize', 12)
  grid on
  fprintf('Axial Sress Plot: Complete\n')
end

% -------------- Shear stress --------------------------------------------------
if strcmp(plotWhat.shearStr, 'yes')
  shearStr = figure('Visible', 'on');
%   figure(4)
  try
    tauSubSet = tau(1:plotWhat.interval:end); % select tau of interest to plot
  catch
    % warning('Failed to limit tau to descrete intervals. Plotting all tau.')
    tauSubSet = tau;
  end

  hold on
  for k = 1:length(tauSubSet)
    plot(rArr*1000, tauSubSet{k}*force*10^3, [marker(k),'-'], 'MarkerIndices', 1:5:length(rArr), 'LineWidth', 1.5)
  end
  plot(rArr*1000, tau{end}*force*10^3, [marker(k+1),'-'], 'MarkerIndices', 1:5:length(rArr), 'LineWidth', 1.5)
  
  grid on
  xlabel(['Radius [', lng_unit, ']'])
  ylabel(['Shear Stress [kPa]'])
  legend(legTxt, 'Location', 'northeast')
%   legend('\omega_{min}, 0 sec', '\omega_{min}, 14 sec', '\omega_{min}, 28 sec', '\omega_{min}, 42 sec', '\omega_{min}, 52 sec',...
% '\omega_{max}, 0 sec', '\omega_{max}, 14 sec', '\omega_{max}, 28 sec', '\omega_{max}, 42 sec', '\omega_{max}, 52 sec', 'NumColumns',2)
  set(gca, 'FontSize', 12)
  fprintf('Shear Stress Plot: Complete\n')
end


% -------------- Peak stress ---------------------------------------------------
if strcmp(plotWhat.peakStr, 'yes')
  peakStr = figure('Visible','on');
  hold on
  
  try
    pl_sub_set = results.peakloc(1:plotWhat.interval:end);
    ps_sub_set = results.peakstr(1:plotWhat.interval:end);
  catch
    pl_sub_set = results.peakloc;
    ps_sub_set = results.peakstr;
  end
  
  yyaxis left; plot(results.vel,pl_sub_set(1)*lng, '-', 'Color', [0 0.4470 0.7410], 'MarkerIndices', 1:10:results.vel, 'LineWidth', 1.5);
  yyaxis right; plot(results.vel,ps_sub_set(1), '-.o', 'MarkerIndices', 1:10:results.vel, 'LineWidth', 1.5);
  for k = 1:length(pl_sub_set)
    yyaxis left; plot(results.vel,pl_sub_set(k)*lng, '-', 'Color', [0 0.4470 0.7410], 'MarkerIndices', 1:10:results.vel, 'LineWidth', 1.5);
    yyaxis right; plot(results.vel,ps_sub_set(k), '-.o', 'MarkerIndices', 1:10:results.vel, 'LineWidth', 1.5);
  end
%   yyaxis right; plot(results.vel, ones(length(results.time)), 'k--', 'LineWidth', 1.5)

  xlabel('Angular velocity [rpm]')
  yyaxis left
  ylabel(['Peak SR Location [', lng_unit, ']'])
  yyaxis right
  ylabel('Strength Ratio')
  legend('Peak SR Location', 'SR value', 'Location', 'southeast')
  set(gca, 'FontSize', 12)
end

% -------------- Max Stress FC  ------------------------------------------------
if strcmp(plotWhat.maxStr, 'yes')
    radMax(1:30) = (results.sArr{1}(3,1:30,1)/mat.stren{1}(3));
    radMax(31:60) = (results.sArr{1}(3,31:60,1)/mat.stren{2}(3));
    radMax(31:60) = (results.sArr{1}(3,31:60,1)/mat.stren{2}(3));
    
    hoopMax(1:30) = (results.sArr{1}(1,1:30,1)/mat.stren{1}(1));
    hoopMax(31:60) = (results.sArr{1}(1,31:end,1)/mat.stren{2}(1));
    
    figure(), hold on
    plot(rArr*1000, radMax)
    plot(rArr*1000, hoopMax)
end

% ------------- Strength Ratio -------------------------------------------------
if strcmp(plotWhat.sr,'yes')
  figure()
  hold on
  
  try
    subSet = results.SR(1:plotWhat.interval:end);
  catch
    subSet = results.SR;
  end

  hold on

  for k = 1:length(subSet)
    plot(rArr*lng, subSet{k}, [marker(k),'-'], 'MarkerIndices', 1:5:length(rArr),  'LineWidth', 1.5);
  end
  plot(rArr*lng, results.SR{end}, [marker(k),'-'], 'MarkerIndices', 1:5:length(rArr),  'LineWidth', 1.5);
  ylabel('Strength Ratio')
  xlabel(['Radius [', lng_unit, ']'])
  legend(legTxt, 'Location', 'northeast')
  grid on
  set(gca, 'Fontsize', 12)

end

%% -----------------------------------------------------------------------------
% Make .gifs
% ------------------------------------------------------------------------------

% -------------- Radial displacement -------------------------------------------
if strcmp(plotWhat.disGif,'yes')
  prog = waitbar(0,'0','Name','Radial Displacement .gif', 'CreateCancelBtn',...
    'setappdata(gcbf,''Canceling'',1)');
  setappdata(prog,'Canceling',0);

  disFig = figure('Visible', 'off');

  [r,t] = meshgrid(rArr, linspace(0,2*pi,length(rArr)));
  x = r .* cos(t);
  y = r .* sin(t);
  z = zeros(length(x),length(y));
  c = ones(length(x),length(y)) .* uArr(1,:) * 1000;

  surf(x,y,z,c, 'EdgeColor', 'none');
  colormap('jet');
  c = colorbar('Eastoutside', 'FontSize', 12);
  c.Label.String = 'Displacement [mm]';

  viscircles(oc, or, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2);

  view(0,90)
  xlabel('Length [m]')
  ylabel('Length [m]')

  currT = 0;
  str = ['Time: ', num2str(currT)];
  loc = [0.613, .621, .3,.3];
  annotation('textbox',loc,'String',str, 'FitBoxToText', 'on');
  set(disFig, 'nextplot', 'replacechildren')

  F = getframe(disFig);
  [frames(:,:,1,1),map] = rgb2ind(F.cdata,256,'nodither');
  for b = 1:vari
    c = ones(length(x),length(y)) .* uArr(b,:) * 1000;

    surf(x,y,z,c, 'EdgeColor', 'none');
    c = colorbar('Eastoutside', 'FontSize', 12);
    c.Label.String = 'Displacement [mm]';

    if mod(b,plotWhat.interval) == 0
      currT = b / plotWhat.interval;
    end
    str = ['Time: ', num2str(currT)];
    annotation('textbox',loc,'String',str, 'FitBoxToText', 'on');

    viscircles(oc, or, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2);

    view(0,90)
    xlabel('Length [m]')
    ylabel('Length [m]')
    set(disFig, 'nextplot', 'replacechildren')

    F = getframe(disFig);
    frames(:,:,1,b) = rgb2ind(F.cdata,map,'nodither');

    perc = (b / vari);
    waitbar(perc,prog,sprintf('%1.0f %%',perc*100))
    if getappdata(prog,'Canceling')
      delete(prog)
      return
    end

  end
  imwrite(frames, map, plotWhat.disGifName, 'DelayTime', plotWhat.delay, 'LoopCount', inf)

  delete(prog)
  fprintf('Radial Displacement gif: Complete\n')
end

% -------------- Radial stress -------------------------------------------------
if strcmp(plotWhat.radGif, 'yes')
  prog = waitbar(0,'0','Name','Radial Stress .gif', 'CreateCancelBtn',...
    'setappdata(gcbf,''Canceling'',1)');
  setappdata(prog,'Canceling',0);

  radFig = figure('Visible','off');

  [r,t] = meshgrid(rArr, linspace(0,2*pi,length(rArr)));
  x = r .* cos(t);
  y = r .* sin(t);
  z = zeros(length(x),length(y));
  c = ones(length(x),length(y)) .* sArr(3,:,1) * 1e-6;

  surf(x,y,z,c, 'EdgeColor', 'none');
  colormap('jet');
  c = colorbar('Eastoutside', 'FontSize', 12);
  c.Label.String = 'Radial Stress [MPa]';

  viscircles(oc, or, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2);

  view(0,90)
  xlabel('Length [m]')
  ylabel('Length [m]')

  currT = 0;
  str = ['Time: ', num2str(currT)];
  loc = [0.613, .621, .3,.3];
  annotation('textbox',loc,'String',str, 'FitBoxToText', 'on');
  set(radFig, 'nextplot', 'replacechildren')

  F = getframe(radFig);
  [frames(:,:,1,1),map] = rgb2ind(F.cdata,256,'nodither');
  for b = 1:vari
    c = ones(length(x),length(y)) .* sArr(3,:,b) * 1e-6;

    surf(x,y,z,c, 'EdgeColor', 'none');
    c = colorbar('Eastoutside', 'FontSize', 12);
    c.Label.String = 'Radial Stress [MPa]';

    if mod(b,plotWhat.interval) == 0
      currT = b / plotWhat.interval;
    end
    str = ['Time: ', num2str(currT)];
    annotation('textbox',loc,'String',str, 'FitBoxToText', 'on');

    viscircles(oc, or, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2);

    view(0,90)
    xlabel('Length [m]')
    ylabel('Length [m]')
    set(radFig, 'nextplot', 'replacechildren')
    F = getframe(radFig);
    frames(:,:,1,b) = rgb2ind(F.cdata,map,'nodither');

    perc = (b / vari);
    waitbar(perc,prog,sprintf('%1.0f %%',perc*100))
    if getappdata(prog,'Canceling')
      delete(prog)
      return
    end
  end
  imwrite(frames, map, plotWhat.radialGifName, 'DelayTime', plotWhat.delay, 'LoopCount', inf)

  delete(prog)
  fprintf('Radial Sress gif: Complete\n')
end

% -------------- Hoop stress ---------------------------------------------------
if strcmp(plotWhat.hoopGif, 'yes')
  prog = waitbar(0,'0','Name','Hoop Stress .gif', 'CreateCancelBtn',...
    'setappdata(gcbf,''Canceling'',1)');
  setappdata(prog,'Canceling',0);
  hoopFig = figure('Visible','off');

  [r,t] = meshgrid(rArr, linspace(0,2*pi,length(rArr)));
  x = r .* cos(t);
  y = r .* sin(t);
  z = zeros(length(x),length(y));
  c = ones(length(x),length(y)) .* sArr(1,:,1) * 1e-6;

  surf(x,y,z,c, 'EdgeColor', 'none');
  colormap('jet');
  c = colorbar('Eastoutside', 'FontSize', 12);
  c.Label.String = 'Hoop Stress [MPa]';

  viscircles(oc, or, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2);

  view(0,90)
  xlabel('Length [m]')
  ylabel('Length [m]')

  currT = 0;
  str = ['Time: ', num2str(currT)];
  loc = [0.613, .621, .3,.3];
  annotation('textbox',loc,'String',str, 'FitBoxToText', 'on');
  set(hoopFig, 'nextplot', 'replacechildren')

  F = getframe(hoopFig);
  [frames(:,:,1,1),map] = rgb2ind(F.cdata,256,'nodither');
  for b = 1:vari
    c = ones(length(x),length(y)) .* sArr(1,:,b) * 1e-6;

    surf(x,y,z,c, 'EdgeColor', 'none');
    c = colorbar('Eastoutside', 'FontSize', 12);
    c.Label.String = 'Hoop Stress [MPa]';

    if mod(b,plotWhat.interval) == 0
      currT = b / plotWhat.interval;
    end
    str = ['Time: ', num2str(currT)];
    annotation('textbox',loc,'String',str, 'FitBoxToText', 'on');

    viscircles(oc, or, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2);

    view(0,90)
    xlabel('Length [m]')
    ylabel('Length [m]')
    set(hoopFig, 'nextplot', 'replacechildren')
    F = getframe(hoopFig);
    frames(:,:,1,b) = rgb2ind(F.cdata,map,'nodither');

    perc = (b / vari);
    waitbar(perc,prog,sprintf('%1.0f %%',perc*100))
    if getappdata(prog,'Canceling')
      delete(prog)
      return
    end
  end
  imwrite(frames, map, plotWhat.hoopGifName, 'DelayTime', plotWhat.delay, 'LoopCount', inf)

  delete(prog)
  fprintf('Hoop Stress gif: Complete\n')
end
