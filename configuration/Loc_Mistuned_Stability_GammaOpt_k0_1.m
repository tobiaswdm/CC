simulation = 'LocalizationSingleSectorStability';

% System parameters
sys.kappa_c = 0.006;
sys.Gamma_Scale = 0.33;

sys.sigma_omega = 0.01;

exc.k = 1;

% Maximum of clearance normalized amplitude
simsetup.LocalizationSingleSectorStability.xi_max = 20;
% Number of samples of clearance normalized amplitude
simsetup.LocalizationSingleSectorStability.Nxi = 2000;
% Range of excitaiton frequencies
simsetup.LocalizationSingleSectorStability.r_range = [0.95 1.02];
% Number of samples of excitation frequencies
simsetup.LocalizationSingleSectorStability.Nr = 4000;
% Number of Mistuning realizations for frequency amplitude curves
simsetup.LocalizationSingleSectorStability.N_MCS = 200;
% Take every stepsize-th point of contour for stability analysis
simsetup.LocalizationSingleSectorStability.stepsize = 100;