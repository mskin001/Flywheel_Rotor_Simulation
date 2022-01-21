function plotStressStrain(rdiv, legTxt, unit)
%% -----------------------------------------------------------------------------
% Define global variables, arrays, and structures
% ------------------------------------------------------------------------------
global rim rArr plotWhat results
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
  rad_data = csvread('C2-1 rad.csv');
  p1 = plot(rad_data(:,1), rad_data(:,2), 'kd', 'MarkerFaceColor', 'k');
  hoop_data = csvread('C2-1 cir.csv');
  p2 = plot(hoop_data(:,1), hoop_data(:,2), 'kd', 'MarkerFaceColor', 'k');
  axi_data = csvread('C2-1 axi.csv');
  p3 = plot(axi_data(:,1), axi_data(:,2), 'kd', 'MarkerFaceColor', 'k');
  
  p4 = plot(rArr(1:30)/rArr(1), sArr{1}(3,1:30)*force/31,'-', 'Color', [0 0.4470 0.7410], 'Linewidth', 1.5);
  p5 = plot(rArr(1:30)/rArr(1), sArr{1}(1,1:30)*force/1062, '-', 'Color', [0.6350 0.0780 0.1840], 'Linewidth', 1.5);
  p6 = plot(rArr(1:30)/rArr(1), sArr{1}(2,1:30)*force/31, '-', 'Color', [0.4660 0.6740 0.1880], 'Linewidth', 1.5);
  p7 = plot(rArr(31:end)/rArr(1), sArr{1}(3,31:end)*force/56,'-', 'Color', [0 0.4470 0.7410], 'Linewidth', 1.5);
  p8 = plot(rArr(31:end)/rArr(1), sArr{1}(1,31:end)*force/3500, '-', 'Color', [0.6350 0.0780 0.1840],  'Linewidth', 1.5);
  p9 = plot(rArr(31:end)/rArr(1), sArr{1}(2,31:end)*force/56, '-', 'Color', [0.4660 0.6740 0.1880], 'Linewidth', 1.5);
  grid on
  legend([p1, p2, p3, p4, p5, p6], {'Ha 1999', 'Ha 1999', 'Ha 1999', 'Model Radial', 'Model Circ.', 'Model Axial'}...
      ,'Location', 'southeast')
  xlabel('Normalized Radius [r/r_1]')
  ylabel('Normalized Stress')
  set(gca, 'FontSize', 12)
  fprintf('Custom plot 1: Complete\n')
end

if strcmp(plotWhat.custom2, 'yes')
    for k = 1:length(results.time)
        inner(1:5,k) = [sArr{k}(1:3,1); tau{k}(1); results.SR{k}(1)];
        outer(1:5,k) = [sArr{k}(1:3,end); tau{k}(end); results.SR{k}(end)];
        rim_inter(1:5,k) = [sArr{k}(1:3,rdiv+1); tau{k}(rdiv+1); results.SR{k}(rdiv+1)];
        hub_inter(1:5,k) = [sArr{k}(1:3,rdiv); tau{k}(rdiv); results.SR{k}(rdiv)];
        rim_mid(1:5,k) = [sArr{k}(1:3, (rdiv+rdiv/2)); tau{k}((rdiv+rdiv/2)); results.SR{k}((rdiv+rdiv/2))];
    end
    % Plots stress at various points against velocity
    figure(), hold on
    plot(results.vel*10^-3, inner(1,:)*10^-6, [marker(1),'-'], 'MarkerIndices', 1:8:length(results.time), 'LineWidth', 1.5)
    plot(results.vel*10^-3, outer(1,:)*10^-6, [marker(2),'-'], 'MarkerIndices', 1:8:length(results.time), 'LineWidth', 1.5)
    plot(results.vel*10^-3, rim_inter(1,:)*10^-6, [marker(3),'-'], 'MarkerIndices', 1:8:length(results.time), 'LineWidth', 1.5)
    plot(results.vel*10^-3, hub_inter(1,:)*10^-6, [marker(4),'-'], 'MarkerIndices', 1:8:length(results.time), 'LineWidth', 1.5)
    plot(results.vel*10^-3, rim_mid(1,:)*10^-6, [marker(5),'-'], 'MarkerIndices', 1:8:length(results.time), 'LineWidth', 1.5)
    xlabel('Angular Velocity [krpm]'), ylabel('Circumferential Stress [MPa]')
    legend('Inner Radius', 'Outer Radius', 'Rim Interface', 'Hub Interface',...
        'Rim Peak Str', 'Location', 'northwest')
    set(gca, 'FontSize', 20), grid on
    figure(), hold on
    plot(results.vel*10^-3, inner(3,:)*10^-6, [marker(1),'-'], 'MarkerIndices', 1:8:length(results.time), 'LineWidth', 1.5)
    plot(results.vel*10^-3, outer(3,:)*10^-6, [marker(2),'-'], 'MarkerIndices', 1:8:length(results.time), 'LineWidth', 1.5)
    plot(results.vel*10^-3, rim_inter(3,:)*10^-6, [marker(3),'-'], 'MarkerIndices', 1:8:length(results.time), 'LineWidth', 1.5)
    plot(results.vel*10^-3, hub_inter(3,:)*10^-6, [marker(4),'-'], 'MarkerIndices', 1:8:length(results.time), 'LineWidth', 1.5)
    plot(results.vel*10^-3, rim_mid(3,:)*10^-6, [marker(5),'-'], 'MarkerIndices', 1:8:length(results.time), 'LineWidth', 1.5)
    xlabel('Angular Velocity [krpm]'), ylabel('Radial Stress [MPa]')
    legend('Inner Radius', 'Outer Radius', 'Rim Interface', 'Hub Interface',...
        'Rim Peak Str', 'Location', 'northwest')
    set(gca, 'FontSize', 12), grid on
    figure(), hold on
    plot(results.vel*10^-3, inner(4,:)*10^-3, [marker(1),'-'], 'MarkerIndices', 1:8:length(results.time), 'LineWidth', 1.5)
    plot(results.vel*10^-3, outer(4,:)*10^-3, [marker(2),'-'], 'MarkerIndices', 1:8:length(results.time), 'LineWidth', 1.5)
    plot(results.vel*10^-3, rim_inter(4,:)*10^-3, [marker(3),'-'], 'MarkerIndices', 1:8:length(results.time), 'LineWidth', 1.5)
    plot(results.vel*10^-3, hub_inter(4,:)*10^-3, [marker(4),'-'], 'MarkerIndices', 1:8:length(results.time), 'LineWidth', 1.5)
    plot(results.vel*10^-3, rim_mid(4,:)*10^-3, [marker(5),'-'], 'MarkerIndices', 1:8:length(results.time), 'LineWidth', 1.5)
    xlabel('Angular Velocity [krpm]'), ylabel('Shear Stress [kPa]')
    legend('Inner Radius', 'Outer Radius', 'Rim Interface', 'Hub Interface',...
        'Rim Peak Str')
    set(gca, 'FontSize', 12), grid on
    figure(), hold on
    plot(results.vel*10^-3, inner(5,:), [marker(1),'-'], 'MarkerIndices', 1:8:length(results.time), 'LineWidth', 1.5)
    plot(results.vel*10^-3, outer(5,:), [marker(2),'-'], 'MarkerIndices', 1:8:length(results.time), 'LineWidth', 1.5)
    plot(results.vel*10^-3, rim_inter(5,:), [marker(3),'-'], 'MarkerIndices', 1:8:length(results.time), 'LineWidth', 1.5)
    plot(results.vel*10^-3, hub_inter(5,:), [marker(4),'-'], 'MarkerIndices', 1:8:length(results.time), 'LineWidth', 1.5)
    plot(results.vel*10^-3, rim_mid(5,:), [marker(5),'-'], 'MarkerIndices', 1:8:length(results.time), 'LineWidth', 1.5)
    xlabel('Angular Velocity [krpm]'), ylabel('SR')
    legend('Inner Radius', 'Outer Radius', 'Rim Interface', 'Hub Interface',...
        'Rim Peak Str')
    set(gca, 'FontSize', 12), grid on
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
  figure(1)
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
%   plot(rArr*1000, tau{end}*force*10^3, [marker(k+1),'-'], 'MarkerIndices', 1:5:length(rArr), 'LineWidth', 1.5)

  grid on
  xlabel(['Radius [', lng_unit, ']'])
  ylabel(['Shear Stress [kPa]'])
  legend(legTxt, 'Location', 'northeast')
  legend('\omega_{min}, 0 sec', '\omega_{min}, 14 sec', '\omega_{min}, 28 sec', '\omega_{min}, 42 sec', '\omega_{min}, 51.5 sec',...
'\omega_{max}, 0 sec', '\omega_{max}, 14 sec', '\omega_{max}, 28 sec', '\omega_{max}, 42 sec', '\omega_{max}, 51.5 sec', 'NumColumns',2)
  set(gca, 'FontSize', 12)
  fprintf('Shear Stress Plot: Complete\n')
end


% -------------- Peak stress ---------------------------------------------------
if strcmp(plotWhat.peakStr, 'yes')
    peakStr = figure('Visible','on');
    hold on
    yyaxis left; plot(results.vel,results.peak_loc*lng, '-', 'Color', [0 0.4470 0.7410], 'MarkerIndices', 1:10:results.vel, 'LineWidth', 1.5);
    yyaxis right; plot(results.vel,results.peak_str, '-.o', 'MarkerIndices', 1:10:results.vel, 'LineWidth', 1.5);

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
    warning('Max Stress FC still under construction')

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
