simulation = 'DetermineSlowFlow';

% System parameters
sys.kappa_c = 0.01;
%sys.kappa_c = 0.2;
sys.Gamma_Scale = 0.15;

exc.k = 1;
H=10;

exc.harmonic.r = 0.99*sqrt(1+4*sys.kappa_c*sin(exc.k*pi/10)^2)/sqrt(0.98);
exc.harmonic.r = 1.0006;
exc.harmonic.r = 1.00008*sqrt(1+4*sys.kappa_c*sin(exc.k*pi/10)^2);

exc.k = 1;
H=10;