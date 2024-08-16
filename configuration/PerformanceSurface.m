% This configuration file starts a time simulation of the performance
% surface and generates the results shown in Fig. 15a-c in the 
% paper T. Weidemann, L. A. Bergman, A. F. Vakakis, M. Krack. (2024)
% "Energy Transfer and Localization in a Forced Cyclic Chain of
% Oscillators with Vibro-Impact Nonlinear Energy Sinks".
% In order to obtain all curves in Fig. 15b and c, the excitation
% wavenumber needs to be adjusted successively

simulation = 'VariationCouplingAndClearance';

%% System parameters
% Tuned system
sys.N_s = 10;               % Number of sectors
sys.kappa_c = 3;            % Linear coupling strength
sys.epsilon_a = 0.02;       % Mass ratio of VI-NES
sys.D = 1e-3;               % Uniform Modal Damping Ratio
sys.Gamma_Scale = 0.4;      % Clearance normalized by linear resonance amplitude
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
exc.k = 0;                  % Excitation wavenumber

% Parameters of nominal coupling
% Range
simsetup.VariationCouplingAndClearance.Range_kappa_c = [5e-3,3];
% Scaling of discrete points in interval
simsetup.VariationCouplingAndClearance.Scaling_kappa_c = 'logarithmic';
% Number of discrete points in interval
simsetup.VariationCouplingAndClearance.Number_kappa_c = 50;

% Parameters of nominal clearance normalized by tuned resonance amplitude
% Range
simsetup.VariationCouplingAndClearance.Range_GammaScale = [0.1, 2];
% Scaling of discrete points in interval
simsetup.VariationCouplingAndClearance.Scaling_GammaScale = 'logarithmic';
% Number of discrete points in interval
simsetup.VariationCouplingAndClearance.Number_GammaScale = 50;

% Frequency stepping range
simsetup.VariationCouplingAndClearance.r_scale = [0.98 1.03];
% Steps in frequency stepping interval
simsetup.VariationCouplingAndClearance.N_rSteps = 50;