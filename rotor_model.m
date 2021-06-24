function [complete] = rotor_model(rim, delta, sigb, mats, compFunc, h, rdiv, profile)

%% -----------------------------------------------------------------------------
% Define global variables
% ------------------------------------------------------------------------------
% Global variables and arrays
global U w arraySize rArr uArr sArr eArr tauArr vari t
% Global structures
global mat results

addpath('ComplianceFunctions')
initial_acc = 0; % rad/s^2
%% -----------------------------------------------------------------------------
% Program Begin
% ------------------------------------------------------------------------------
b = 1;
[~, cols] = size(profile);

while b <= cols
    w = (pi/30) * profile(2,b); %initial angular velocity
%     w0 = w;

    t = profile(1,b);
    %   fprintf('Create Variable Arrays: Complete\n')
    %% ---------------------------------------------------------------------------
    % Preallocate variables
    % ----------------------------------------------------------------------------
    arraySize = length(rim);
    U = zeros(1,arraySize);
    rArr = zeros(1,(arraySize-1)*rdiv);    % radius vector for descretization
    uArr = zeros(1,(arraySize-1)*rdiv);    % displacement vector for discretization
    sArr = zeros(4,(arraySize-1)*rdiv);    % stress vector
    tauArr = zeros(1,(arraySize-1)*rdiv);
    eArr = zeros(4, rdiv);    % strain vector in each direction
    E = zeros(1,cols);
    %   fprintf('Preallocate Memory: Complete\n')

    %% ---------------------------------------------------------------------------
    % Create Q matrices for all materials
    % ----------------------------------------------------------------------------
    for k = 1:length(mats)
        func = compFunc{k};
        mat.file{k} = ['MaterialProperties\', mats{k}];
        matProp = load(mat.file{k});
        mat.Q{1,k} = stiffMat(matProp.mstiff, func);
        mat.rho{k} = matProp.rho;

        try
          mat.stren{k} = matProp.stren;
        catch
          break
        end
    end

    %   fprintf('Create Material Property Matrices: Complete\n')
    %% ---------------------------------------------------------------------------
    % Calculate displacement magnitude at the inner and outer surface of each rim
    % these are used as boundary conditions to find C. ~ is used to disregard
    % output of force vector results. These can be important for debugging and
    % verification purposes, but are not necessary for the program. Check function
    % discription for mor info
    [~, ~, ~, ~] = boundaryConditions(sigb, delta);

    %   fprintf('Calculate Boundary Conditions: Complete\n')
    %% ---------------------------------------------------------------------------
    % Calculate discrete displacement, stain, and stress for each rim ~ here is
    % used to the [C] matrix output. This is useful for debugging and
    % verification purposes but not necessary for the function. Check function
    % description for mor info
    [~] = discretizeStressStrain(rdiv, delta);

    %   fprintf('Descretize Stress/Strain: Complete\n')

    %% ----------------------------------------------------------------------------
    % Calculate the share stress on the rim.
    if b == 1
        alpha = initial_acc;
    else
        alpha = (profile(2,b-1) - profile(2,b)) / (profile(1,b-1) - profile(1,b));
    end

    [~] = shearStress(alpha, rdiv);
    
    %% ---------------------------------------------------------------------------
    % Calculate the current stored energy in the rotor based on the current
    % angluar velocity
    
    E(b) = find_energy(h);
    %% ---------------------------------------------------------------------------
    % Store results for post processing
    % ----------------------------------------------------------------------------
    results.uArr{b} = uArr;
    results.sArr{b} = sArr;
    results.tauArr{b} = tauArr;
    results.vel(b) = w * (30 / pi);
    %   fprintf('Current time: %5.2f\n', b*tStep)
    %   fprintf('Iteration %2.0f Complete\n', b)

    b = b + 1;

end

%% -----------------------------------------------------------------------------
% Calculate failure criterion
% ------------------------------------------------------------------------------
vari = cols;
[SR] = failureIndex(rdiv);
results.SR = SR;
results.E = E;
complete = 1;
% 
% % fprintf('Create Output Plots: Complete\n\n')
% fprintf('Program Complete\n')
