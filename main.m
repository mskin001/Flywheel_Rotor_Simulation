clc
clear
close('all','force')
format long
%% -----------------------------------------------------------------------------
% Define global variables
% ------------------------------------------------------------------------------
% Global variables and arrays
global rim
% Global structures
global plotWhat results

%% -----------------------------------------------------------------------------
% Define initial conditions and rotor size
% ------------------------------------------------------------------------------
% Rotor
r1 = .17518; %linspace(6.875, 10, 5) * 0.0254;
r2 = .17463; %linspace(6.897, 10, 5) * 0.0254;
r3 = .3175; %linspace(12.5, 15, 5) * 0.0254;

% rim = [0.0762, 0.09144, 0.10668]; %Tzeng2001 
% rim = [r1, r2, r3];
rdiv = 30; % number of points per rim to analyze
% delta = [.000378, 0]; % Tzeng 2012 press fit [m]
% delta = [.0004, .0004, 0]; % m
delta = [.000254, 0]; % m

sigb = [0, 0]; % [Pa]
% mats = {'IM7_8552_Tzeng2001.mat', 'IM7_8552_Tzeng2001.mat'};
mats = {'Walkingshaw_GFRP_withFoS.mat' 'Walkingshaw_CFRP_withFoS.mat'};
compFunc = {'no', 'no'}; % compliance function, input 'no' to turn off creep modeling
h = .12573; % rotor thickness in [m]

% Speed/velocity
profile = [1;...           % [ t1 t2 t3;
           3125];             %   v1 v2 v3]

% Plotting
% legTxt = {'Current model', 'Aparicio 2011'};
legTxt = {'0 sec', '1 year', '5 years'}; % Controls legend entries for graphs
plotWhat.custom1 = 'no';        % any custom plot. Go to plotStressStrain.m to modify (first if statement)
plotWhat.maxStr = 'yes';        % maximum stress failure criteria
plotWhat.radDis = 'no';          % Radial displacement v. radius
plotWhat.radStr = 'yes';         % Radial stress v. radius plot
plotWhat.hoopStr = 'yes';        % Hoop stress v. radius plot
plotWhat.shearStr = 'no';       % Shear stress v. radius
plotWhat.peakStr = 'no';        % 2-yaxis plot. Peak stress location and SR v. time
plotWhat.sr = 'yes';

plotWhat.disGif = 'no';          % Displacement gif, surface plot
plotWhat.disGifName = 'Displacement.gif';
plotWhat.radGif = 'no';          % Radial stress gif, surface plot
plotWhat.radialGifName = 'Radial Stress.gif';
plotWhat.hoopGif = 'no';         % Hoop stress gif, surface plot
plotWhat.hoopGifName = 'Hoop Stress.gif';

plotWhat.interval = 1;          % Display time interval on figures
plotWhat.delay = 0;             % Time delay in seconds between frames in the gifs,
                                % 0 is fastest
                                
E = zeros(length(r1), length(r2), length(r3));
sr = zeros(length(r1), length(r2), length(r3));
%% -----------------------------------------------------------------------------
% Start Program
% ------------------------------------------------------------------------------
fprintf('Number of rims: %2.0f\n',length(rim)-1)
fprintf('Material Selections: %s\n', mats{1:end})
fprintf('Max rotational velocity: %6.3f rpm\n\n', profile(end,end))
fprintf('Program Start: \n')

for k = 1:length(r1)
    rim(1) = r1(k);
   for j = 1:length(r2)
        rim(2) = rim(1) + r2(j);
      for x = 1:length(r3)
            rim(3) = rim(1) + r3(x);
            complete = rotor_model(rim, delta, sigb, mats, compFunc, h, rdiv, profile);
            sr(k,j,x) = results.SR(1);
            E(k,j,x) = results.E(1);
            disp(k), disp(j), disp(x)
      end
      
   end
%     file_name = ['exp', num2str(1)];
%     rad_disp(k,:) = results.uArr{1};
%     rad_str(k,:) = results.sArr{1}(3,:,1);
%     hoop_str(k,:) = results.sArr{1}(1,:,1);
    
    
end
%% -----------------------------------------------------------------------------
% Make Plots
% ------------------------------------------------------------------------------
plotStressStrain(legTxt)