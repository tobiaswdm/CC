%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Cylic chain of oscialltors with vibro-impact absorbers
%
% Code written by Tobias Weidemann, M.Sc. - (C) 2024
% University of Stuttgart, Germany
% Institute of Aircraft Propulsion Systems
%
% Contact: tobias.weidemann@ila.uni-stuttgart.de
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clearvars;

% Randomize Seed
rng("shuffle");

%% Configuration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% System configuration file in configuration folder
configuration = 'VariationCouplingAndClearanceMCS_eigenfreqmist_k0_1';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add system paths
fprintf('Adding folder paths... \n')
addpath('.\configuration\') % Path for config files
addpath('.\style\')         % Path for style functions
DefaultStyle;
addpath('.\system\')        % Path for system functions
addpath('.\simulations\')   % Path for simulation cases
addpath('.\integrator\')    % Path for Time integration schemes
addpath('.\postprocessing\')% Path for Post processings
addpath('.\statistics\')    % Path for Statistics
addpath('.\analytics\')     % Path for Analytical methods
addpath('.\harmonic_balance\')     % Path for HB
fprintf('Initializing data structures... \n')
InitStructs;

% Create path to save data and figures
savepath = ['.\data\' configuration '\'];
savepath_backup = ['.\data\' configuration '\backup\'];
[~,~] = mkdir(savepath);
[~,~] = mkdir(savepath_backup);

% Load configuration file
fprintf('Loading configuration... \n')
run([configuration '.m'])

% Run desired simulation
switch simulation
    case 'VariationCouplingAndClearanceMCS'
        fprintf(['Running Variation of Coupling and Clearance' ...
            ' with MCS... \n'])
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
        fprintf(['Running analytical study on localization in a single' ...
            ' sector... \n'])
        LocalizationSingleSectorAnalytical;
    case 'LocalizationSingleSectorStability'
        fprintf(['Running stability analysis on localization in a' ...
            ' single sector... \n'])
        LocalizationSingleSectorStability;
    case 'GsaprStability'
        fprintf('Running stability analysis of GSAPR... \n')
        GsaprStability;
    case 'DetermineSlowFlow'
        fprintf(['Solving Slow Flow equation for AQPR in frequency' ...
            ' domain... \n'])
        DetermineSlowFlow;
    case 'SlowFlowTimeSimulation'
        fprintf('Solving Slow Flow equation for AQPR in time domain... \n')
        SlowFlowTimeSimulation;
    case 'LinearMistuningAnalysis'
        fprintf('Performing Linear Mistuning Analysis using MCS... \n')
        LinearMistuningAnalysis;
    otherwise
        fprintf('Oops... How did we end up here? \n')
        error('Simulation not defined')
end

[~,~] = rmdir(savepath_backup,'s');