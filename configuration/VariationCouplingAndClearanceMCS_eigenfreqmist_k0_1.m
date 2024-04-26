simulation = 'VariationCouplingAndClearanceMCS';

exc.k = 1;
sys.N_s = 10;

sys.sigma_omega = 1e-2;

% Number of Monte Carlos Simulations per nominal configuration
simsetup.VariationCouplingAndClearanceMCS.N_MCS = 22000;

% Parameters of nominal coupling
% Range
simsetup.VariationCouplingAndClearanceMCS.Range_kappa_c = [1e-3,1];
% Scaling of discrete points in interval
simsetup.VariationCouplingAndClearanceMCS.Scaling_kappa_c = 'logarithmic';
% Number of discrete points in interval
simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c = 4;

% Parameters of nominal clearance normalized by tuned resonance amplitude
% Range
simsetup.VariationCouplingAndClearanceMCS.Range_GammaScale = [0.2,0.55];
% Scaling of discrete points in interval
simsetup.VariationCouplingAndClearanceMCS.Scaling_GammaScale = 'logarithmic';
% Number of discrete points in interval
simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale = 3;

% Frequency stepping range
simsetup.VariationCouplingAndClearanceMCS.r_scale = [0.98 1.03];
% Steps in frequency stepping interval
simsetup.VariationCouplingAndClearanceMCS.N_rSteps = 50;