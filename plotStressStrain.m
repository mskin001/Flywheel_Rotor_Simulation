function plotStressStrain(legTxt)
%% -----------------------------------------------------------------------------
% Define global variables, arrays, and structures
% ------------------------------------------------------------------------------
global rim rArr plotWhat results mat
uArr = results.uArr;
sArr = results.sArr;
tau =  results.tauArr;

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
  plot(rArr*39.3701, sArr{1}(3,:,1)*.000145038,'-', 'Linewidth', 1.5)
  rad_data = csvread('radial-stress_t1_pressfit_cylinder.csv');
  plot(rad_data(:,1), rad_data(:,2), 'kd', 'MarkerFaceColor', 'k')
  
  plot(rArr*39.3701, sArr{2}(3,:,1)*.000145038,'--', 'Linewidth', 1.5)
  rad_data = csvread('radial-stress_t2_pressfit_cylinder.csv');
  plot(rad_data(:,1), rad_data(:,2), 'ks', 'MarkerFaceColor', 'k')
  
  plot(rArr*39.3701, sArr{3}(3,:,1)*.000145038,':', 'Linewidth', 1.5)
  rad_data = csvread('radial-stress_t3_pressfit_cylinder.csv');
  plot(rad_data(:,1), rad_data(:,2), 'k^', 'MarkerFaceColor', 'k')
  
  grid on
  xlabel('Radius [in]')
  ylabel('Stress [psi]')
  legend('Model t1', 'Tzeng t1', 'Model t10', 'Tzeng t10', 'Model tinf', 'Tzeng tinf'...
      ,'Location', 'southeast')
  set(gca, 'FontSize', 12)
  
  figure(), hold on
  plot(rArr*39.3701, sArr{1}(1,:,1)*.000145038, '-', 'Linewidth', 1.5)
  hoop_data = csvread('hoop-stress_t1_pressfit_cylinder.csv');
  plot(hoop_data(:,1), hoop_data(:,2), 'kd', 'MarkerFaceColor', 'k')
    
  plot(rArr*39.3701, sArr{2}(1,:,1)*.000145038, '--', 'Linewidth', 1.5)
  hoop_data = csvread('hoop-stress_t2_pressfit_cylinder.csv');
  plot(hoop_data(:,1), hoop_data(:,2), 'ks', 'MarkerFaceColor', 'k')
  
  plot(rArr*39.3701, sArr{3}(1,:,1)*.000145038, ':', 'Linewidth', 1.5)
  hoop_data = csvread('hoop-stress_t3_pressfit_cylinder.csv');
  plot(hoop_data(:,1), hoop_data(:,2), 'k^', 'MarkerFaceColor', 'k')
  
%   plot(rArr*1000, tau{1}*10^-6, ':', 'Color', [0.4940 0.1840 0.5560], 'Linewidth', 1.5)
%   plot(tau_data(:,1)*1000, tau_data(:,2), 'ko')
  grid on
  legend('Model t1', 'Tzeng t1', 'Model t10', 'Tzeng t10', 'Model tinf', 'Tzeng tinf'...
      ,'Location', 'southeast')
  xlabel('Radius [in]')
  ylabel('Stress [psi]')
%   legend('Model 10*\sigma_r', 'Aparicio2011, 10*\sigma_r','Model \sigma_\theta',...
%       'Aparicio2011, \sigma_\theta', 'Model \tau_r_\theta', 'Aparicio2011, \tau_r_\theta',...
%       'NumColumns', 3, 'Location', 'southoutside')
  set(gca, 'FontSize', 12)
  fprintf('Custom plot 1: Complete\n')
  
%   figure(), hold on
%   plot(rArr*39.3701, uArr{1}*39.3701, 'Linewidth', 1.5)
%   disp_data = csvread('Tzeng2001RadialDispInitial.csv');
%   plot(disp_data(:,1), disp_data(:,2), 'kd')
%   
%   plot(rArr*39.3701, uArr{2}*39.3701, 'Linewidth', 1.5)
%   disp_data = csvread('Tzeng2001RadialDispMid.csv');
%   plot(disp_data(:,1), disp_data(:,2), 'ks')
%   
%   plot(rArr*39.3701, uArr{3}*39.3701, 'Linewidth', 1.5)
%   disp_data = csvread('Tzeng2001RadialDispInf.csv');
%   plot(disp_data(:,1), disp_data(:,2), 'k^')
%   grid on
%   xlabel('Radius [in]')
%   ylabel('Stress [psi]')
%   legend('Model t1', 'Tzeng t1', 'Model t10', 'Tzeng t10', 'Model tinf', 'Tzeng tinf')
%   set(gca, 'FontSize', 12)
end

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
%% -----------------------------------------------------------------------------
% Make figures
% ------------------------------------------------------------------------------

% -------------- Radial displacements ------------------------------------------
if strcmp(plotWhat.radDis, 'yes')
  radDis = figure('Visible','on');
  hold on
  try
    subSet = uArr(plotWhat.interval:plotWhat.interval:end);
  catch
    subSet = uArr;
  end

  hold on
  for k = 1:length(subSet)
    plot(rArr*39.3701, subSet{k}(1,:)*39.3701, 'LineWidth', 1.5);
  end
  xlabel('Radius [in]')
  ylabel('Radial Displacement [in]')
  legend(legTxt, 'Location', 'southeast')
%   legend('Tzeng Initial', 'Tzeng 10 years', 'Tzeng Infinite', 'Initial','10 Years', 'Infinite')
  set(gca, 'FontSize', 12)
  grid on
  fprintf('Radial Displacement Plot: Complete\n')
end

% -------------- Radial stress -------------------------------------------------
if strcmp(plotWhat.radStr, 'yes')
  radStr = figure('Visible','on');
  hold on


  try
    subSet = sArr(plotWhat.interval:plotWhat.interval:end);
  catch
    subSet = sArr;
  end

  hold on
%   plot(rArr*1000, sArr{1}(3,:,1)*10^-6, 'LineWidth', 1.5)
  for k = 1:length(subSet)
    plot(rArr*39.3701, subSet{k}(3,:,1)*10^-6, 'LineWidth', 1.5);
  end
  xlabel('Radius [in]')
  ylabel('Radial Stress [MPa]')
  legend(legTxt, 'Location', 'southeast')
%   legend('Tzeng Initial', 'Tzeng 10 years', 'Tzeng Infinite', 'Initial','10 Years', 'Infinite')
  set(gca, 'FontSize', 12)
  grid on
  fprintf('Radial Sress Plot: Complete\n')
end

% -------------- Hoop stress ---------------------------------------------------
if strcmp(plotWhat.hoopStr, 'yes')
  hoopStr = figure('Visible','on'); %#ok<*NASGU>
  hold on

  try
    subSet = sArr(plotWhat.interval:plotWhat.interval:end);
  catch
    subSet = sArr;
  end

  hold on
%   plot(rArr*1000, sArr{1}(1,:,1)*10^-6, 'LineWidth', 1.5)
  for k = 1:length(subSet)
    plot(rArr*39.3701, subSet{k}(1,:,1)*10^-6, 'LineWidth', 1.5);
  end


  xlabel('Radius [in]')
  ylabel('Circumferential Stress [MPa]')
  legend(legTxt, 'Location', 'southeast', 'NumColumns', 5)
  set(gca, 'FontSize', 12)
  grid on
  fprintf('Hoop Stress Plot: Complete\n')
end

% -------------- Shear stress --------------------------------------------------
if strcmp(plotWhat.shearStr, 'yes')
  shearStr = figure('Visible', 'on');
  try
    tauSubSet = tau(plotWhat.interval:plotWhat.interval:end); % select tau of interest to plot
  catch
    % warning('Failed to limit tau to descrete intervals. Plotting all tau.')
    tauSubSet = tau;
  end

  hold on
  for k = 1:1
    plot(rArr*1000, tau{end}, 'LineWidth', 1.5)
  end
%   stressData = csvread('aparicio2011_results.csv', 1, 0);
%   plot(stressData(:,1)*1000, stressData(:,2), 'k*')
% %   for k = 1:length(tauSubSet)
%     plot(rArr*1000,tauSubSet{k}, 'LineWidth', 1.5);
%   end
  xlabel('Radius [mm]')
  ylabel('Shear Stress [Pa]')
%   legend(legTxt, 'Location', 'northeast')
  set(gca, 'FontSize', 12)
  fprintf('Shear Stress Plot: Complete\n')
end


% -------------- Peak stress ---------------------------------------------------
if strcmp(plotWhat.peakStr, 'yes')
  peakStr = figure('Visible','on');
  hold on
  yyaxis left; plot(results.vel,results.peakloc*1000, '-', 'Color', [0 0.4470 0.7410], 'MarkerIndices', 1:10:results.vel, 'LineWidth', 1.5);
  yyaxis right; plot(results.vel,results.peakstr, '-.o', 'MarkerIndices', 1:10:results.vel, 'LineWidth', 1.5);
%   yyaxis right; plot(results.vel, ones(length(results.time)), 'k--', 'LineWidth', 1.5)

  xlabel('Angular velocity [rpm]')
  yyaxis left
  ylabel('Peak SR Location [mm]')
  yyaxis right
  ylabel('Strength Ratio')
  legend('Peak SR Location', 'SR value', 'Location', 'southeast')
  set(gca, 'FontSize', 12)
end

% ------------- Strength Ratio --------------------------------------------
if strcmp(plotWhat.sr,'yes')
  figure()
  hold on
  
  plot(rArr*39.3701, results.SR{1}, '-o', 'Color', [0 0.4470 0.7410],...
      'MarkerIndices', 1:5:length(rArr), 'Linewidth', 1.5)
  plot(rArr*39.3701, results.SR{2}, '--d','Color', [0.6350 0.0780 0.1840],...
      'MarkerIndices', 1:5:length(rArr), 'Linewidth', 1.5)
  plot(rArr*39.3701, results.SR{3}, ':v', 'Color', [0.4940 0.1840 0.5560],...
      'MarkerIndices', 1:5:length(rArr), 'Linewidth', 1.5)
  
  ylabel('Strength Ratio')
  xlabel('Radius [in]')
%   legend('SR \it t=1', 'SR \it t=10', 'SR \it t=20')
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
