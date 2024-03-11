simulation = 'ConvergenceAMCS';

%% System parameters
sys.kappa_c = 0.01;
sys.Gamma_Scale = 0.33;    % Gamma = qref*Gamma_scale

% Mistuned system
sys.sigma_omega = 0.01;
sys.sigma_g = 0;

exc.k = 1;

% Number of Monte Carlos Simulations per nominal configuration
% Maximum Number of MCSs
simsetup.ConvergenceAMCS.N_MCS = 100;
% Number of MCS Permutations for converegence analysis
simsetup.ConvergenceAMCS.Shuffle_MCS = 30;
% Minimum Number of MCSs in evaluation of AMCS
simsetup.ConvergenceAMCS.Min_MCS = 30;


% Frequency stepping range
simsetup.ConvergenceAMCS.r_scale = [0.98 1.03];
% Steps in frequency stepping interval
simsetup.ConvergenceAMCS.N_rSteps = 50;