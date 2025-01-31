% This configuration file starts a time simulation of a SSMR and generates
% the results shown in Fig. 21-middle in the paper
% T. Weidemann, L. A. Bergman, A. F. Vakakis, M. Krack. (2025)
% "Energy transfer and localization in a forced cyclic chain of
% oscillators with vibro-impact nonlinear energy sinks".
% doi: https://doi.org/10.1007/s11071-025-10928-4

simulation = 'FindCoexistingSolutions';

%% System parameters
% Tuned system
sys.N_s = 10;               % Number of sectors
sys.kappa_c = 0.002;        % Linear coupling strength
sys.epsilon_a = 0.02;       % Mass ratio of VI-NES
sys.D = 1e-3;               % Uniform Modal Damping Ratio
sys.Gamma_Scale = 0.33;     % Clearance normalized by linear resonance amplitude
sys.eN = 0.8;               % Restitution coefficient

%% Excitation
exc.type = 'harmonic';      % Excitation type 'harmonic' or 'sweep'
exc.k = 1;                  % Excitation wavenumber
exc.harmonic.r = 0.992297;

% Number of Randomized Initial Conditions that are tested
simsetup.FindCoexisitingSolutions.N_MCS = 1e4;