simulation = 'SynchronizationSingleSectorAnalytical';

% System parameters
sys.kappa_c = 0.006;
sys.Gamma_Scale = 0.33;

sys.sigma_omega = 0.01;

exc.k = 1;

% Maximum of clearance normalized amplitude
simsetup.SynchronizationSingleSectorAnalytical.xi_max = 50;
% Number of samples of clearance normalized amplitude
simsetup.SynchronizationSingleSectorAnalytical.Nxi = 3000;
% Range of excitaiton frequencies
simsetup.SynchronizationSingleSectorAnalytical.r_range = [0.93 1.07];
% Number of samples of excitation frequencies
simsetup.SynchronizationSingleSectorAnalytical.Nr = 3000;
% Number of Mistuning realizations for frequency amplitude curves
simsetup.SynchronizationSingleSectorAnalytical.N_MCS = 200;