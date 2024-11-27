% This configuration file tracks the bifurcation point of the first
% resonance branch of the LSR with a synchroniztation in one sector and
% generates the results shown in Fig. 18 in the paper
% T. Weidemann, L. A. Bergman, A. F. Vakakis, M. Krack. (2024)
% "Energy Transfer and Localization in a Forced Cyclic Chain of
% Oscillators with Vibro-Impact Nonlinear Energy Sinks".
% For all wavenumbers, exc.k needs to be adjusted accordingly.

simulation = 'FollowBifurcationFRSofLSR';

%% System parameters
% Tuned system
sys.N_s = 10;               % Number of sectors
sys.epsilon_a = 0.02;       % Mass ratio of VI-NES
sys.D = 1e-3;               % Uniform Modal Damping Ratio
sys.eN = 0.8;               % Restitution coefficient
sys.Gamma_Scale = 0.33;     % Normalized clearance

%% Excitation
exc.k = 0;                  % Excitation wavenumber

% Parameters of clearance-normalized amplitude
% Maximum of clearance normalized amplitude
simsetup.FollowBifurcationFRSofLSR.xi_max = 20;
% Number of samples of clearance normalized amplitude
simsetup.FollowBifurcationFRSofLSR.Nxi = 4000;

% Parameters of excitation frequency
% Sampling range
simsetup.FollowBifurcationFRSofLSR.r_range = [0.9 1.1];
% Number of samples
simsetup.FollowBifurcationFRSofLSR.Nr = 3000;

% Parameters of nominal coupling
% Range
simsetup.FollowBifurcationFRSofLSR.kappac_range = [1e-3,1];
% Number of discrete points in interval
simsetup.FollowBifurcationFRSofLSR.Nkappac = 50;
% Step Size for stability of frequency amplitude curve
simsetup.FollowBifurcationFRSofLSR.stepsize = 1;