function [residual] = SlowFlowResidual(coeff_ReIm,H,k_mod,sys,exc)
%SLOWFLOWRESIDUAL Get the Harmonic Balance residual for the 
% almost quasi periodic response (AQPR)
%
% H - Harmonic order
% coeff_ReIm - [2*(2*H+1),1] Real and Imaginary parts of complex Fourier
%                coefficients of slow flow
%               [Phi_-H_Re; ...; Phi_0_Re; ...; Phi_H_Re;
%                Phi_-H_Im; ...; Phi_0_Im; ...; Phi_H_Im]
%
% H - Harmonic order
% k_mod - Wavenumber of slow modulation

% Fast excitation frequency
r = exc.harmonic.r;

% Modulation wavenumber phase angle
theta_k_mod = 2*pi*k_mod/sys.N_s;

% Excitation wavenumber phase angle
theta_k_0 = 2*pi*exc.k/sys.N_s;

% Convert to complex Fourier Coefficients
Phi_hat = [coeff_ReIm(1:(H+1));10;coeff_ReIm((H+2):(2*H))]...
+1i*coeff_ReIm((2*H)+(1:(2*H+1)));


% Transform to time domain
% chose at least 1000 sampling points per slow period
N = max(1000,4*(H+1));

% Evaluate nonlinear dispersion relation for average amplitude level
% get slow modulation frequency Omega (group velocity)
%[~,Omega_d] = NonlinearDispersionRelation(abs(Phi_hat(H+1)),theta_k_0,sys);
Omega = coeff_ReIm(end)/1000;

% Evaluate the real an imaginary part of the slow nonlinear force
% in the time domain and get Fourier coefficients of dynamic components 
% force coefficients without scaling of -8*epsilon_a*r^2 / pi^2
[F_hat] = EvaluateSlowForceCoefficient(Phi_hat,N,H,sys);


% Build residual equation

% Harmonics
h = (-H:H)';

% Sectors
j = 0:(sys.N_s-1);

% Complex stiffness of Harmonics h
Sh = -(1-sys.epsilon_a)*r^2 + 1 + 2*sys.kappa_c * ...
     (1-cos(theta_k_mod*h + theta_k_0)) ....
     - 2*r*Omega*(1-sys.epsilon_a)*h ...
     + 1i*r*(exp(1i*(theta_k_mod*h + theta_k_0)*j)*sys.C(:,1));

% Complex residual
complex_residual = Sh.*Phi_hat - ...
                   8*sys.epsilon_a*(r/pi)^2 * F_hat - ...
                   circshift(eye(2*H+1,1),H,1);



% Split into real and imaginary component
residual = [real(complex_residual); ...
            imag(complex_residual)];

%residual(H+2) = coeff_ReIm(H+2)-20;

end