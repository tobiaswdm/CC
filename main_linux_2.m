clc;
close all;
clearvars;


%% Add system paths
fprintf('Adding folder paths... \n')
addpath('./configuration/') % Path for config files
addpath('./style/')         % Path for style functions
DefaultStyle;
addpath('./system/')        % Path for system functions
addpath('./simulations/')   % Path for simulation cases
addpath('./integrator/')    % Path for Time integration schemes
addpath('./postprocessing/')% Path for Post processings
addpath('./statistics/')    % Path for Statistics
addpath('./analytics/')     % Path for Analytical methods
addpath('./harmonic_balance/')     % Path for HB
fprintf('Initializing data structures... \n')
InitStructs;

% Randomize Seed
rng("shuffle");

%% Configuration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

configuration = 'ConvMCS_kc0_01_sigmag_0_01_sigmaomega_0_01_Gamma_0_33_k0_1';  % System configuration file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create path to save data and figures
savepath = [savepath configuration '/'];
[~,~] = mkdir(savepath);

savepath = ['/data1/tweidemann/CCOO_mt/data/' configuration '/'];
savepath_backup = ['/data1/tweidemann/CCOO_mt/data/' configuration '/backup/'];
[~,~] = mkdir(savepath);
[~,~] = mkdir(savepath_backup);

% Load configuration file
fprintf('Loading configuration... \n')
run([configuration '.m'])

% Run desired simulation
switch simulation
    case 'VariationCouplingAndClearanceMCS'
        fprintf('Running Variation of Coupling and Clearance with MCS... \n')
        VariationCouplingAndClearanceMCS;
    case 'VariationCouplingAndClearance'
        fprintf('Running Variation of Coupling and Clearance... \n')
        VariationCouplingAndClearance;
    case 'ConvergenceMCS'
        fprintf('Running convergence study on MCS... \n')
        ConvergenceMCS;
    case 'TimeSimulation'
        fprintf('Running single time simulation... \n')
        TimeSimulation;
    case 'LocalizationSingleSectorAnalytical'
        fprintf('Running analytical study on localization in a single sector... \n')
        LocalizationSingleSectorAnalytical;
    case 'LocalizationSingleSectorStability'
        fprintf('Running stability analysis on localization in a single sector... \n')
        LocalizationSingleSectorStability;
    case 'DetermineSlowFlow'
        fprintf('Solving Slow Flow equation for AQPR... \n')
        DetermineSlowFlow;
    otherwise
        fprintf('Oops... How did we end up here? \n')
        error('Simulation not defined')
end


[~,~] = rmdir(savepath_backup,'s');