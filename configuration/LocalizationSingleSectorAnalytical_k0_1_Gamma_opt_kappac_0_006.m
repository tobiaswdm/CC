simulation = 'LocalizationSingleSectorAnalytical';

% System parameters
sys.kappa_c = 0.006;
sys.Gamma_Scale = 0.33;

sys.sigma_omega = 0.01;

exc.k = 1;

% Maximum of clearance normalized amplitude
simsetup.LocalizationSingleSectorAnalytical.xi_max = 50;
% Number of samples of clearance normalized amplitude
simsetup.LocalizationSingleSectorAnalytical.Nxi = 3000;
% Range of excitaiton frequencies
simsetup.LocalizationSingleSectorAnalytical.r_range = [0.93 1.07];
% Number of samples of excitation frequencies
simsetup.LocalizationSingleSectorAnalytical.Nr = 3000;
% Number of Mistuning realizations for frequency amplitude curves
simsetup.LocalizationSingleSectorAnalytical.N_MCS = 200;