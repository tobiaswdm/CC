% This configuration file starts a stability analysis 
% for a LSR with synchronization in a single sector and generates the 
% results shown in Fig. 3, Fig. 10, Fig. 11 and Fig. 13 in the paper
% T. Weidemann, L. A. Bergman, A. F. Vakakis, M. Krack. (2025)
% "Energy transfer and localization in a forced cyclic chain of
% oscillators with vibro-impact nonlinear energy sinks".
% doi: https://doi.org/10.1007/s11071-025-10928-4

simulation = 'SynchronizationSingleSectorStability';

%% System parameters
% Tuned system
sys.N_s = 10;               % Number of sectors
sys.kappa_c = 0.05;         % Linear coupling strength
sys.epsilon_a = 0.02;       % Mass ratio of VI-NES
sys.D = 1e-3;               % Uniform Modal Damping Ratio
sys.Gamma_Scale = 0.07;     % Clearance normalized by linear resonance amplitude
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
simsetup.SynchronizationSingleSectorStability.xi_max = 20;
% Number of samples of clearance normalized amplitude
simsetup.SynchronizationSingleSectorStability.Nxi = 6000;
% Range of excitaiton frequencies
simsetup.SynchronizationSingleSectorStability.r_range = [0.75 1.1];
% Number of samples of excitation frequencies
simsetup.SynchronizationSingleSectorStability.Nr = 4000;
% Number of Mistuning realizations for frequency amplitude curves
simsetup.SynchronizationSingleSectorStability.N_MCS = 0;
% Take every stepsize-th point of contour for stability analysis
simsetup.SynchronizationSingleSectorStability.stepsize = 150;