% This configuration file plots the FRS of the LSR with synchronization in
% a single sector and generates the results shown in Fig. 12a of the paper 
% T. Weidemann, L. A. Bergman, A. F. Vakakis, M. Krack. (2025)
% "Energy transfer and localization in a forced cyclic chain of
% oscillators with vibro-impact nonlinear energy sinks".
% doi: https://doi.org/10.1007/s11071-025-10928-4

simulation = 'SynchronizationSingleSectorAnalytical';

%% System parameters
% Tuned system
sys.N_s = 10;               % Number of sectors
sys.kappa_c = 0.05;         % Linear coupling strength
sys.epsilon_a = 0.02;       % Mass ratio of VI-NES
sys.D = 1e-3;               % Uniform Modal Damping Ratio
sys.Gamma_Scale = 0.15;     % Clearance normalized by linear resonance amplitude
sys.eN = 0.8;               % Restitution coefficient

% Mistuned system
sys.sigma_omega = 0;        % Rel. Standard deviation of local eigefrequencies
sys.sigma_g = 0;            % Rel. Standard deviation of local clearances
sys.adjustC = false;        % Set true if sys.D should also refer to
                            % modes of the mistuned system and false if
                            % damping matrix of tuned system should also be
                            % used for the mistuned case


%% Excitation
exc.type = 'harmonic';      % Excitation type 'harmonic' or 'sweep'
exc.k = 1;                  % Excitation wavenumber

%% Simulation parameters
% Maximum of clearance normalized amplitude
simsetup.SynchronizationSingleSectorAnalytical.xi_max = 50;
% Number of samples of clearance normalized amplitude
simsetup.SynchronizationSingleSectorAnalytical.Nxi = 4000;
% Range of excitaiton frequencies
simsetup.SynchronizationSingleSectorAnalytical.r_range = [0.93 1.12];
% Number of samples of excitation frequencies
simsetup.SynchronizationSingleSectorAnalytical.Nr = 2000;
% Number of Mistuning realizations for frequency amplitude curves
simsetup.SynchronizationSingleSectorAnalytical.N_MCS = 0;