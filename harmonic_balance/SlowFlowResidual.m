function [residual] = SlowFlowResidual(Phi_hat_ReIm,H,k_mod,sys,exc)
%SLOWFLOWRESIDUAL Get the Harmonic Balance residual for the 
% almost quasi periodic response (AQPR)
%
% Phi_hat_ReIm - [2*H+1,1] Real and Imaginary parts of complex Fourier
%                coefficients of nonlinear force
%
% H - Hamronic order
% k_mod  -Wavenumber of slow modulation


% Modulation wavenumber phase angle
theta_k_mod = 2*pi*k_mod/sys.N_s;

% Excitation wavenumber phase angle
theta_k_0 = 2*pi*exc.k/sys.N_s;

% Convert to complex Fourier Coefficients
Phi_hat = [Phi_hat_ReIm(1);...
    Phi_hat_ReIm(1+(1:H))+1i*Phi_hat_ReIm((1+H)+(1:H))];

% Transform to time domain
% chose at least 1000 sampling points per slow period
N = max(1000,4*(H+1));
% Slow flow in time domain
slow_flow = FrequencyTime(Phi_hat,N,'Freq_to_Time');

% Evaluate the real an imaginary part of the slow nonlinear force
% in the time domain
[force_real,force_imag]=...
    EvaluateSlowForceCoefficient(slow_flow,sys);

% Get Fourier coefficients of dynamic components force coefficients
% without scaling of 8*epsilon_a*r^2 / pi^2

% Fundamental of real component
F_hat_real_0 = mean(force_real);

% All harmonics of imaginary component
F_hat_imag = FrequencyTime(force_imag,H,'Time_to_Freq');

% Evaluate nonlinear dispersion relation and group velocity to get
% fast excitaiton frequency r and the slow modulation frequency Omega
[r,Omega] = ...
    NonlinearDispersionRelation(Phi_hat_ReIm(1),F_hat_real_0,sys,exc);

% Build residual equation

% Harmonics
h = (0:H)';

% Complex residual
complex_stiffness = sys.C(1,1)+2*1i*(r*Omega*(1-sys.epsilon_a)*h-...
    sys.kappa_c*sin(theta_k_mod*h)*sin(theta_k_0));

% Complex residual
complex_residual = complex_stiffness.*Phi_hat +...
        8*sys.epsilon_a*r^2/pi^2 * F_hat_imag - eye(H+1,1);

% Split into real and imaginary component
residual = [real(complex_residual(1:(H+1))); ...
            imag(complex_residual(2:(H+1)))];

end