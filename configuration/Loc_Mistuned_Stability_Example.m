simulation = 'SynchronizationSingleSectorStability';

% System parameters
sys.kappa_c = 0.05;
sys.Gamma_Scale = 0.08;

sys.sigma_omega = 0.01;

exc.k = 1;

% Maximum of clearance normalized amplitude
simsetup.SynchronizationSingleSectorStability.xi_max = 20;
% Number of samples of clearance normalized amplitude
simsetup.SynchronizationSingleSectorStability.Nxi = 4000;
% Range of excitaiton frequencies
simsetup.SynchronizationSingleSectorStability.r_range = [0.8 1.1];
% Number of samples of excitation frequencies
simsetup.SynchronizationSingleSectorStability.Nr = 4000;
% Number of Mistuning realizations for frequency amplitude curves
simsetup.SynchronizationSingleSectorStability.N_MCS = 0;
% Take every stepsize-th point of contour for stability analysis
simsetup.SynchronizationSingleSectorStability.stepsize = 150;