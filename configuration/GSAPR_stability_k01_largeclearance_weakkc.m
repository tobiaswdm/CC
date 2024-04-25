%simulation = 'DetermineSlowFlow';
simulation = 'GsaprStability';
% System parameters
sys.kappa_c = 0.01;
sys.Gamma_Scale = 0.3;

exc.k = 1;

% Maximum of clearance normalized amplitude
simsetup.GsaprStability.xi_max = 20;
% Number of samples of clearance normalized amplitude
simsetup.GsaprStability.Nxi = 1000;
% Range of excitaiton frequencies
simsetup.GsaprStability.r_range = [0.96 1.03]*...
    sqrt(1+4*sys.kappa_c*sin(exc.k*pi/10)^2);
% Number of samples of excitation frequencies
simsetup.GsaprStability.Nr = 1000;
% Take every stepsize-th point of contour for stability analysis
simsetup.GsaprStability.stepsize = 60;

