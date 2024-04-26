simulation = 'ConvergenceMCS';

%% System parameters
sys.kappa_c = 0.01;
sys.Gamma_Scale = 0.33;    % Gamma = qref*Gamma_scale

% Mistuned system
sys.sigma_omega = 0.01;
sys.sigma_g = 0.01;

exc.k = 2;

% Number of Monte Carlos Simulations per nominal configuration
% Maximum Number of MCSs
simsetup.ConvergenceMCS.N_MCS = 100000;
% Number of MCS Permutations for converegence analysis
simsetup.ConvergenceMCS.Shuffle_MCS = 2000;

% Frequency stepping range
simsetup.ConvergenceMCS.r_scale = [0.98 1.03];
% Steps in frequency stepping interval
simsetup.ConvergenceMCS.N_rSteps = 50;