% This configuration file starts a time simulation of a stability analysis 
% for a LSR with synchronization in a single sector and generates parts of 
% the results shown in Fig. 19 of the paper
% T. Weidemann, L. A. Bergman, A. F. Vakakis, M. Krack. (2024)
% "Energy Transfer and Localization in a Forced Cyclic Chain of
% Oscillators with Vibro-Impact Nonlinear Energy Sinks".

simulation = 'SynchronizationSingleSectorStability';

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
exc.k = 0;                  % Excitation wavenumber

%% Simulation Setup
% Maximum of clearance normalized amplitude
simsetup.SynchronizationSingleSectorStability.xi_max = 10;
% Number of samples of clearance normalized amplitude
simsetup.SynchronizationSingleSectorStability.Nxi = 4000;
% Range of excitaiton frequencies
simsetup.SynchronizationSingleSectorStability.r_range = [0.8 1.1];
% Number of samples of excitation frequencies
simsetup.SynchronizationSingleSectorStability.Nr = 2000;
% Number of Mistuning realizations for frequency amplitude curves
simsetup.SynchronizationSingleSectorStability.N_MCS = 0;
% Take every stepsize-th point of contour for stability analysis
simsetup.SynchronizationSingleSectorStability.stepsize = 95;