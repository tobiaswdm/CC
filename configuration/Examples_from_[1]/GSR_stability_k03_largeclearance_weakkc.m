% This configuration file starts a stability analysis of the GSR for a
% clearance that is close to the bifurcation point of the FRS and generates
% the results shown in Fig. 8-bottom-right of the paper 
% T. Weidemann, L. A. Bergman, A. F. Vakakis, M. Krack. (2025)
% "Energy transfer and localization in a forced cyclic chain of
% oscillators with vibro-impact nonlinear energy sinks".
% doi: https://doi.org/10.1007/s11071-025-10928-4

simulation = 'GSRStability';

%% System parameters
% Tuned system
sys.N_s = 10;               % Number of sectors
sys.kappa_c = 0.01;         % Linear coupling strength
sys.epsilon_a = 0.02;       % Mass ratio of VI-NES
sys.D = 1e-3;               % Uniform Modal Damping Ratio
sys.Gamma_Scale = 0.3;      % Clearance normalized by linear resonance amplitude
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
exc.k = 3;                  % Excitation wavenumber

%% Simulation Setup
% Maximum of clearance normalized amplitude
simsetup.GSRStability.xi_max = 20;
% Number of samples of clearance normalized amplitude
simsetup.GSRStability.Nxi = 1000;
% Range of excitaiton frequencies
simsetup.GSRStability.r_range = [0.9 1.03]*...
    sqrt(1+4*sys.kappa_c*sin(exc.k*pi/10)^2);
% Number of samples of excitation frequencies
simsetup.GSRStability.Nr = 1000;
% Take every stepsize-th point of contour for stability analysis
simsetup.GSRStability.stepsize = 60;