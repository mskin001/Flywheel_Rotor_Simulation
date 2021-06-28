function plotStressStrain(legTxt)
%% -----------------------------------------------------------------------------
% Define global variables, arrays, and structures
% ------------------------------------------------------------------------------
global rim rArr plotWhat results mat
uArr = results.uArr;
sArr = results.sArr;
tau =  results.tauArr;

face = ['s', 'd', '^', 'o', '*', 'x', '<', '>', 'v', '+', '|', 'p', 'h'];
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
  
    for k = 1:length(results.sArr)
        radMax(k) = max(sArr{k}(3,:,1))/mat.stren{1}(3);
        hoopMax(k) = max(sArr{k}(1,:,1))/mat.stren{1}(1);
        tauMax(k) = max(tau{k}) / mat.stren{1}(end);
    end
        
    hold on
    plot(results.vel, radMax,'-', 'Color', [0 0.4470 0.7410], 'Linewidth', 1.5)
    %   plot(rad_data(:,1), rad_data(:,2), 'kv')
    plot(results.vel, hoopMax, '--', 'Color', [0.6350 0.0780 0.1840], 'Linewidth', 1.5)
    %   plot(hoop_data(:,1), hoop_data(:,2), 'k^')
%       plot(results.vel, tauMax, ':', 'Color', [0.4940 0.1840 0.5560], 'Linewidth', 1.5)
    %   plot(tau_data(:,1)*1000, tau_data(:,2), 'ko')

    xlabel('Angular Velocity [rpm]')
    ylabel('Normalized Stress')
    legend('Radial', 'Circ.', 'Location', 'northwest')
    %   legend('Model 10*\sigma_r', 'Aparicio2011, 10*\sigma_r','Model \sigma_\theta',...
    %       'Aparicio2011, \sigma_\theta', 'Model \tau_r_\theta', 'Aparicio2011, \tau_r_\theta',...
    %       'NumColumns', 3, 'Location', 'southoutside')
    grid on, set(gca, 'FontSize', 12)

    figure(), hold on
    plot(results.vel, radMax,'-', 'Color', [0 0.4470 0.7410], 'Linewidth', 1.5)
    plot(results.vel, hoopMax, '--', 'Color', [0.6350 0.0780 0.1840], 'Linewidth', 1.5)
    plot(results.vel, results.peakstr, ':s', 'Color', [0.4940 0.1840 0.5560],...
        'MarkerIndices', 1:3:length(results.vel), 'MarkerFaceColor', [0.4940 0.1840 0.5560], 'Linewidth', 1.5)
    xlabel('Angular Velocity [rpm]')
    ylabel('Normalized Stress')
    legend('Radial', 'Circ.', 'Peak SR', 'Location', 'northwest')
    grid on,  set(gca, 'FontSize', 12)
    fprintf('Custom plot 1: Complete\n')
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
    plot(rArr*1000, subSet{k}(1,:), 'LineWidth', 1.5);
  end
  xlabel('Radius [mm]')
  ylabel('Radial Displacement [m]')
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
  plot(rArr*1000, sArr{1}(3,:,1)*10^-6, 'Marker', face(1), 'MarkerIndices', 1:4:length(rArr), 'LineWidth', 1.5)
  for k = 1:length(subSet)
    plot(rArr*1000, subSet{k}(3,:,1)*10^-6, 'Marker', face(k+1), 'MarkerIndices', 1:4:length(rArr), 'LineWidth', 1.5);
  end
  plot(rArr*1000, sArr{end}(3,:,1)*10^-6, 'Marker', face(k+2), 'MarkerIndices', 1:4:length(rArr), 'LineWidth', 1.5)
  
  xlabel('Radius [mm]')
  ylabel('Radial Stress [MPa]')
  legend(legTxt, 'Location', 'southeast', 'NumColumns', 2)
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
  plot(rArr*1000, sArr{1}(1,:,1)*10^-6, 'Marker', face(1), 'MarkerIndices', 1:4:length(rArr), 'LineWidth', 1.5)
  for k = 1:length(subSet)
    plot(rArr*1000, subSet{k}(1,:,1)*10^-6, 'Marker', face(k+1), 'MarkerIndices', 1:4:length(rArr), 'LineWidth', 1.5);
  end
  plot(rArr*1000, sArr{end}(1,:,1)*10^-6, 'Marker', face(k+2), 'MarkerIndices', 1:4:length(rArr), 'LineWidth', 1.5)

  xlabel('Radius [mm]')
  ylabel('Circumferential Stress [MPa]')
  legend(legTxt, 'Location', 'northeast', 'NumColumns', 5)
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
  plot(rArr*1000, tau{1}*10^-3, 'Marker', face(1), 'MarkerIndices', 1:4:length(rArr), 'LineWidth', 1.5)
  for k = 1:length(subSet)
    plot(rArr*1000, tauSubSet{k}*10^-3, 'Marker', face(k+1), 'MarkerIndices', 1:4:length(rArr), 'LineWidth', 1.5)
  end
  plot(rArr*1000, tau{end}*10^-3, 'Marker', face(k+2), 'MarkerIndices', 1:4:length(rArr), 'LineWidth', 1.5)

  xlabel('Radius [mm]')
  ylabel('Shear Stress [kPa]')
  legend(legTxt, 'Location', 'northeast')
  grid on, set(gca, 'FontSize', 12)
  fprintf('Shear Stress Plot: Complete\n')
end


% -------------- Peak stress ---------------------------------------------------
if strcmp(plotWhat.peakStr, 'yes')
  peakStr = figure('Visible','on');
  hold on
  yyaxis left; plot(results.vel,results.peakloc*1000, '-', 'Color', [0 0.4470 0.7410], 'MarkerIndices', 1:10:results.vel, 'LineWidth', 1.5);
  yyaxis right; plot(results.vel,results.peakstr, '-.o', 'MarkerIndices', 1:4:results.vel, 'LineWidth', 1.5);
  yyaxis right; plot(results.vel, ones(length(results.time)), 'k--', 'LineWidth', 1.5)

  xlabel('Angular velocity [rpm]')
  yyaxis left
  ylabel('Peak SR Location [mm]')
  yyaxis right
  ylabel('Strength Ratio')
  legend('Peak SR Location', 'SR value', 'Location', 'northwest')
  grid on, set(gca, 'FontSize', 12)
end

% ------------- Strength Ratio --------------------------------------------
if strcmp(plotWhat.sr,'yes')
    figure()
    hold on

    %   plot(rArr*1000, results.SR(1,:), '-o', 'Color', [0 0.4470 0.7410],...
    %       'MarkerIndices', 1:5:length(rArr), 'Linewidth', 1.5)
    %   plot(rArr*1000, results.SR(12,:), '--d','Color', [0.6350 0.0780 0.1840],...
    %       'MarkerIndices', 1:5:length(rArr), 'Linewidth', 1.5)
    %   plot(rArr*1000, results.SR(18,:), '-.^','Color', [0.6350 0.0780 0.1840],...
    %       'MarkerIndices', 1:5:length(rArr), 'Linewidth', 1.5)
    %   plot(rArr*1000, results.SR(end,:), ':v', 'Color', [0.4940 0.1840 0.5560],...
    %       'MarkerIndices', 1:5:length(rArr), 'Linewidth', 1.5)
    try
        subSet = results.SR(plotWhat.interval:plotWhat.interval:end,:);
    catch
        subSet = SR;
    end

    plot(rArr*1000, results.SR(1,:), 'Marker', face(1), 'MarkerIndices', 1:4:length(rArr), 'LineWidth', 1.5)
    for k = 1:min(size(subSet))
        plot(rArr*1000, subSet(k,:), 'Marker', face(k+1), 'MarkerIndices', 1:4:length(rArr), 'LineWidth', 1.5);
    end
    plot(rArr*1000, results.SR(end,:), 'Marker', face(k+2), 'MarkerIndices', 1:4:length(rArr), 'LineWidth', 1.5)

    ylabel('Strength Ratio')
    xlabel('Radius [mm]')
    %   legend('SR \it t=1s', 'SR \it t=2s', 'SR \it t=3s', 'SR \it t=4s', 'SR \it t=5s')
    legend(legTxt, 'Location', 'SouthOutside', 'NumColumns', 5)
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
