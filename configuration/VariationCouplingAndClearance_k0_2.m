simulation = 'VariationCouplingAndClearance';

exc.k = 2;
sys.N_s = 10;

%% System parameters
% Parameters of nominal coupling
% Range
simsetup.VariationCouplingAndClearance.Range_kappa_c = [1e-3,3];
% Scaling of discrete points in interval
simsetup.VariationCouplingAndClearance.Scaling_kappa_c = 'logarithmic';
% Number of discrete points in interval
simsetup.VariationCouplingAndClearance.Number_kappa_c = 50;

% Parameters of nominal clearance normalized by tuned resonance amplitude
% Range
simsetup.VariationCouplingAndClearance.Range_GammaScale = [0.1,1];
% Scaling of discrete points in interval
simsetup.VariationCouplingAndClearance.Scaling_GammaScale = 'logarithmic';
% Number of discrete points in interval
simsetup.VariationCouplingAndClearance.Number_GammaScale = 48;

% Frequency stepping range
simsetup.VariationCouplingAndClearance.r_scale = [0.98 1.03];
% Steps in frequency stepping interval
simsetup.VariationCouplingAndClearance.N_rSteps = 50;