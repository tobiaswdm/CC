% This configuration file starts a time simulation of a stability analysis 
% for a LSR with synchronization in two opposing sectors and generates 
% parts of the results shown in Fig. 20-middle of the paper
% T. Weidemann, L. A. Bergman, A. F. Vakakis, M. Krack. (2025)
% "Energy transfer and localization in a forced cyclic chain of
% oscillators with vibro-impact nonlinear energy sinks".
% doi: https://doi.org/10.1007/s11071-025-10928-4

simulation = 'SynchronizationOpposingSectorsStability';

%% System parameters
% Tuned system
sys.N_s = 10;               % Number of sectors
sys.kappa_c = 0.002;        % Linear coupling strength
sys.epsilon_a = 0.02;       % Mass ratio of VI-NES
sys.D = 1e-3;               % Uniform Modal Damping Ratio
sys.Gamma_Scale = 0.33;     % Clearance normalized by linear resonance amplitude
sys.eN = 0.8;               % Restitution coefficient

% Mistuned system
sys.sigma_omega = 0;        % Rel. Standard deviation of local eigefrequencies
sys.sigma_g = 0;            % Rel. Standard deviation of local clearances
sys.adjustC = false;        % Set true if sys.D should also refer to
                            % modes of the mistuned system and false if
                            % damping matrix of tuned system should also be
                            % used for the mistuned case

%% Excitation
exc.k = 1;                  % Excitation wavenumber

%% Simulation Setup
% Maximum of clearance normalized amplitude
simsetup.SynchronizationOpposingSectorsStability.xi_max = 100;
% Number of samples of clearance normalized amplitude
simsetup.SynchronizationOpposingSectorsStability.Nxi = 6000;
% Range of excitaiton frequencies
simsetup.SynchronizationOpposingSectorsStability.r_range = [0.8 1.1];
% Number of samples of excitation frequencies
simsetup.SynchronizationOpposingSectorsStability.Nr = 2000;
% Take every stepsize-th point of contour for stability analysis
simsetup.SynchronizationOpposingSectorsStability.stepsize = 150;