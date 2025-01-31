% This configuration file starts a time simulation of a ASMR with weak
% inter-sector coupling and generates the results shown in Figs. 6 and 7  
% in the paper
% T. Weidemann, L. A. Bergman, A. F. Vakakis, M. Krack. (2025)
% "Energy transfer and localization in a forced cyclic chain of
% oscillators with vibro-impact nonlinear energy sinks".
% doi: https://doi.org/10.1007/s11071-025-10928-4

simulation = 'BackBoneStability';

%% System parameters
% Tuned system
sys.N_s = 10;               % Number of sectors
sys.epsilon_a = 0.02;       % Mass ratio of VI-NES
sys.D = 1e-3;               % Uniform Modal Damping Ratio
sys.Gamma_Scale = 0.43;     % Clearance normalized by linear resonance amplitude
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

%% Simulation Setup
% Range of inter-sector coupling strengths
simsetup.BackBoneStability.kappac_range = [10^-4, 10^0];
% Number of inter-sector coupling strengths
simsetup.BackBoneStability.Nkappac = 40;
% Number of clearance-normlaized amplitudes
simsetup.BackBoneStability.Nxi = 30;