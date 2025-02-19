%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CC (pronounced Sisi) performs numerical and/or analytical analyses on a 
% Cylic Chain of Oscillators with Vibro-Impact Nonlinear Energy Sinks
%
% The Code for CC was written by:
% Tobias Weidemann - (C) 2024
% University of Stuttgart, Germany
% Institute of Aircraft Propulsion Systems
%
% Contact: tobias.weidemann@ila.uni-stuttgart.de
%
% Feel free to use, share and modify under the GPL-3.0 license.
% If you use CC, please refer to the paper:
%
% T. Weidemann, L. A. Bergman, A. F. Vakakis, M. Krack. (2025)
% "Energy transfer and localization in a forced cyclic chain of
% oscillators with vibro-impact nonlinear energy sinks".
% Nonlinear Dynamics. doi: https://doi.org/10.1007/s11071-025-10928-4
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clearvars;

%% Configuration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% System configuration file in configuration folder
configuration = 'FRS_GSR';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%% ....Starting Simulation with CC.... %%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('-------------------------------------------------------------------')
disp('---------------- Code written by Tobias Weidemann  ----------------')
disp('-- Feel free to use, share and modify under the GPL-3.0 license  --')
disp('-------------------------------------------------------------------')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%% ....Starting Simulation with CC.... %%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')

%% Add system paths
fprintf('Adding folder paths... \n')
addpath(genpath('.\configuration\')) % Path for config files
addpath('.\style\')         % Path for style functions

DefaultStyle;
addpath('.\system\')            % Path for system functions
addpath('.\simulations\')       % Path for simulation cases
addpath('.\integrator\')        % Path for Time integration schemes
addpath('.\postprocessing\')    % Path for Post processings
addpath('.\analytics\')         % Path for Analytical methods
fprintf('Initializing data structures... \n')
InitStructs;

% Randomize Seed
rng("shuffle");

if exist(configuration,'file') 
    % Create path to save data and figures
    savepath = ['.\data\' configuration '\'];
    savepath_backup = ['.\data\' configuration '\backup\'];
    [~,~] = mkdir(savepath);
    [~,~] = mkdir(savepath_backup);
else
    error('Configuration doesn''t exist.')
end

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
    case 'SynchronizationSingleSectorAnalytical'
        fprintf(['Running analytical study on localization in a single' ...
            ' sector... \n'])
        SynchronizationSingleSectorAnalytical;
    case 'SynchronizationSingleSectorStability'
        fprintf(['Running stability analysis on LSR with synchronization in a' ...
            ' single sector... \n'])
        SynchronizationSingleSectorStability;
    case 'GSRStability'
        fprintf('Running stability analysis of GSR... \n')
        GSRStability;
    case 'LinearMistuningAnalysis'
        fprintf('Performing Linear Mistuning Analysis using MCS... \n')
        LinearMistuningAnalysis;
    case 'FollowBifurcationFRSofLSR'
        fprintf('Following the first bifurcation point of the FRS of the LSR... \n')
        FollowBifurcationFRSofLSR;
    case 'SynchronizationOpposingSectorsStability'
        fprintf(['Running stability analysis on LSR with synchronization in ' ...
            'two opposing sectors... \n'])
        SynchronizationOpposingSectorsStability;
    case 'BackBoneStability'
        fprintf('Running stability analysis along Backbone of GSR... \n')
        BackBoneStability;
    case 'SlowFlowTimeSimulation'
        fprintf('Running Slow Flow Simulation... \n')
        SlowFlowTimeSimulation;
    case 'FindCoexistingSolutions'
        fprintf('Finding coexisiting solutions by varying initial conditions... \n')
        FindCoexistingSolutions;
    case 'FrequencyStepping'
        fprintf('Determining Frequency-Amplitude curve through sine-stepping... \n')
        FrequencyStepping;
    otherwise
        fprintf('Oops... How did we end up here? \n')
        error('Simulation not defined')
end

[~,~] = rmdir(savepath_backup,'s');