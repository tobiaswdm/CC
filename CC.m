%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CC (pronounced Sisi) is a tool that performs numerical and/or
% analytical analyses on a Cylic Chain of Oscillators with Vibro-Impact
% Nonlinear Energy Sinks (VI-NESs)
%
% The Code for CC was written by:
% Tobias Weidemann - (C) 2024
% University of Stuttgart, Germany
% Institute of Aircraft Propulsion Systems
%
% Contact: tobias.weidemann@ila.uni-stuttgart.de
%
% Feel free to use, share and modify under the GPL-3.0 license.
% CC is purely academic and comes with no warranty.
% If you use CC for your own research, please refer to the paper:
%
% T. Weidemann, L. A. Bergman, A. F. Vakakis, M. Krack. (2024)
% "Energy Transfer and Localization in a Forced Cyclic Chain of
% Oscillators with Vibro-Impact Nonlinear Energy Sinks".
% Manuscript submitted to Nonlinear Dynamics
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clearvars;

%% Configuration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% System configuration file in configuration folder
configuration = 'VariationCouplingAndClearance_k0_0';

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
addpath('.\configuration\') % Path for config files
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

if isfile(['.\configuration\' configuration '.m']) 
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
    otherwise
        fprintf('Oops... How did we end up here? \n')
        error('Simulation not defined')
end

[~,~] = rmdir(savepath_backup,'s');