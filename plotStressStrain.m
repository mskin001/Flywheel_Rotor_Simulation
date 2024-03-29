function plotStressStrain(results)
%% -----------------------------------------------------------------------------
% Define global variables, arrays, and structures
% ------------------------------------------------------------------------------
global rArr sArr uArr vari
global plotWhat mat

legendLabels = {'\it t = 0, \omega = 0', '\it t = 0', '\it t = 1 yr', '\it t = 5 yrs'};
%% -----------------------------------------------------------------------------
% Define rim origional centers and radii
% ------------------------------------------------------------------------------
numRims = length(plotWhat.rims);
oc = zeros(numRims,2);
or = plotWhat.rims;

%% -----------------------------------------------------------------------------
% Create a custom plot
% ------------------------------------------------------------------------------
% To create additional custom plots copy and paste this section to create
% as many custom plots as desired.

if strcmp(plotWhat.custom1, 'yes')
  haRadData = csvread('Ha99_GFRP_optimized_radialStress.csv');
  hahoopData = csvread('Ha99_GFRP_optimized_hoopStress.csv');
  
  nr = rArr ./ min(rArr);
  
  srRad = sArr(3,:,1) ./ mat.stren{1}(3);
  srHoop = sArr(1,:,1) ./ mat.stren{1}(1);
  
  sw_rad_data = csvread('glass rim-3D Plane Strain-Results-Radial Stress-1.csv',...
    9,1,[9,1,35,2]);
  sw_hoop_data = csvread('glass rim-3D Plane Strain-Results-Hoop Stress-3.csv',...
    9,1,[9,1,35,2]);
  
  sw_norm_rad = [sw_rad_data(:,1)/mat.stren{1}(3),sw_rad_data(:,2)/min(sw_rad_data(:,2))];
  sw_norm_hoop = [sw_hoop_data(:,1)/mat.stren{1}(1),sw_hoop_data(:,2)/min(sw_hoop_data(:,2))];
  
  radStr = figure();
  hold on
  plot(haRadData(:,1),haRadData(:,2),'kv-', 'MarkerFaceColor', 'k')
  plot(hahoopData(:,1),hahoopData(:,2),'k^-', 'MarkerFaceColor', 'k')
  plot(sw_norm_rad(:,2),sw_norm_rad(:,1),'g*', 'MarkerFaceColor', 'g')
  plot(sw_norm_hoop(:,2),sw_norm_hoop(:,1),'gd', 'MarkerFaceColor', 'g')
  plot(nr,srRad,'b--s', 'LineWidth', 1)
  plot(nr,srHoop, 'r--o', 'LineWidth', 1)
  
  axialStr = figure();
  plot(rArr,sArr(2,:,1),'bo-');
%   axis([0.5, 1, -1, 3]);
  grid on
  xlabel('r/r_{min}')
  ylabel('Normalized Stress')
  set(gca, 'FontSize', 12)
  legend('Ha 1999 Radial', 'Ha 1999 Circumfrential','Model Radial','Model Circumfrential', 'Location', 'SouthEast')
  
  fprintf('Custom plot 1: Complete\n')
end

% -----------------------------------------------------------------------------
% Plot displacement
% ------------------------------------------------------------------------------
if strcmp(plotWhat.radDis, 'yes')
  radDis = figure('Visible','on');
  hold on
  
  plot(rArr*1000,uArr(1,:),'b-o', 'LineWidth', 1)
%   plot(rArr*1000,uArr(2,:),'r-d', 'LineWidth', 1)
%   plot(rArr*1000,uArr(3,:),'k-s', 'LineWidth', 1) 
  
  xlabel('Radius [mm]')
  ylabel('Radial Displacement [m]')
  legend(legendLabels, 'Location', 'southeast')
  set(gca, 'FontSize', 12)
  grid on
  fprintf('Radial Displacement Plot: Complete\n')
end

% -----------------------------------------------------------------------------
% Plot radial stress
% ------------------------------------------------------------------------------
if strcmp(plotWhat.radStr, 'yes')
  radStr = figure('Visible','on');
  hold on
  
  try
    zeroVelRadStr = results.sArr(3,:,1)*10^-6;
    plot(rArr*1000,zeroVelRadStr, '--', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2)
  catch
    %do nothing
  end
  plot(rArr*1000,sArr(3,:,1)*10^-6,'b-o', 'LineWidth', 1)
  plot(rArr*1000,sArr(3,:,2)*10^-6,'r-d', 'LineWidth', 1)
  plot(rArr*1000,sArr(3,:,3)*10^-6,'k-s', 'LineWidth', 1)
  
  xlabel('Radius [mm]')
  ylabel('Radial Stress [MPa]')
  legend(legendLabels, 'Location', 'southeast')
  set(gca, 'FontSize', 12)
  grid on
  fprintf('Radial Sress Plot: Complete\n')
end

% -----------------------------------------------------------------------------
% Plot hoop stress
% ------------------------------------------------------------------------------
if strcmp(plotWhat.hoopStr, 'yes')
  hoopStr = figure('Visible','on'); %#ok<*NASGU>
  hold on
  try
    zeroVelHoopStr = results.sArr(1,:,1)*10^-6;  
    plot(rArr*1000,zeroVelHoopStr, '--', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2)
  catch
    %do nothing
  end
  plot(rArr*1000,sArr(1,:,1)*10^-6,'b-o', 'LineWidth', 1)
%   plot(rArr*1000,sArr(1,:,2)*10^-6,'r-d', 'LineWidth', 1)
%   plot(rArr*1000,sArr(1,:,3)*10^-6,'k-s', 'LineWidth', 1)

  xlabel('Radius [mm]')
  ylabel('Circumferential Stress [MPa]')
  legend(legendLabels, 'Location', 'southeast')
  set(gca, 'FontSize', 12)
  grid on
  fprintf('Hoop Stress Plot: Complete\n')
end






%% ------------------------------------------------------------------------
% Make .gifs
% -------------------------------------------------------------------------

% ------ Displacement .gifs -----------------------------------------------
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

% ------ Radial Str .gif --------------------------------------------------
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

% ------ Hoop str .gif ----------------------------------------------------
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










