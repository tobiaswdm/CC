%simulation = 'DetermineSlowFlow';
simulation = 'SlowFlowTimeSimulation';
% System parameters
sys.kappa_c = 0.01;
%sys.kappa_c = 0.2;
sys.Gamma_Scale = 0.15;

exc.k = 1;

sol.N_Tau = 1600;

exc.harmonic.r = sqrt(1+4*sys.kappa_c*sin(exc.k*pi/10)^2);
exc.harmonic.r = 1.00195;

exc.k = 1;
H=10;