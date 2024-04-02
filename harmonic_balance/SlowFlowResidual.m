function [residual] = SlowFlowResidual(coeff_ReIm,H,k_mod,sys,exc)
%SLOWFLOWRESIDUAL Get the Harmonic Balance residual for the 
% almost quasi periodic response (AQPR)
%
% H - Hamronic order
% coeff_ReIm - [2*(2*H+1),1] Real and Imaginary parts of complex Fourier
%                coefficients of nonlinear force
%               [Phi_0; Phi_1_Re; Phi_2_Re; ...; Phi_H_Re;
%                Phi_1_Im; Phi_2_Im; ...; Phi_H_Im;
%                Gamma_0; Gamma_1_Re; Gamma_2_Re; ...; Gamma_H_Re;
%                Gamma_1_Im; Gamma_2_Im; ...; Gamma_H_Im]
%
% H - Hamronic order
% k_mod - Wavenumber of slow modulation


% Modulation wavenumber phase angle
theta_k_mod = 2*pi*k_mod/sys.N_s;

% Excitation wavenumber phase angle
theta_k_0 = 2*pi*exc.k/sys.N_s;

% Convert to complex Fourier Coefficients
Phi_hat = coeff_ReIm(1:(2*H+1))+1i*coeff_ReIm((2*H+1)+(1:(2*H+1)));


% Transform to time domain
% chose at least 1000 sampling points per slow period
N = max(1000,4*(H+1));
% Slow flow in time domain
Phi = FrequencyTime(Phi_hat,N,'Freq_to_Time');

% Evaluate nonlinear dispersion relation and group velocity to get
% slow modulation frequency Omega
[~,Omega] = NonlinearDispersionRelation(coeff_ReIm(1),theta_k_0,sys);

% Evaluate the real an imaginary part of the slow nonlinear force
% in the time domain and get Fourier coefficients of dynamic components 
% force coefficients without scaling of 8*epsilon_a*r^2 / pi^2
[Fnl_hat_real,Fnl_hat_imag]=...
    EvaluateSlowForceCoefficient(Phi_hat,N,H,sys);

% Evaluate the real an imaginary part of the slow coupling forces
% in the time domain and get Fourier coefficients of dynamic components 
% force coefficients without scaling of kappa_c
[Fc_hat_real,Fc_hat_imag]=...
    EvaluateSlowCouplingForceCoefficient(Phi_hat,Gamma_hat,theta_k_mod,...
                                        theta_k_0,N,H);

% Evaluate the real an imaginary part of the forcing term
% in the time domain and get Fourier coefficients
[Ff_hat_real,Ff_hat_imag]=...
    EvaluateSlowExcitationForceCoefficient(Gamma_hat,N,H);

% Evaluate the real an imaginary part of the derivative term without
%sclaing of 2*(1-epsilon_a)
Fd_hat_real= ...
    EvaluateSlowDerivativeForceCoefficient(Phi_hat,Gamma_hat,N,H);


% Build residual equation

% Harmonics
h = (0:H)';

% Linear component stiffnesses
sRe = repmat((-exc.harmonic.r^2*(1-sys.epsilon_a) + 1 + 2*sys.kappa_c),...
    [H+1,1]);
sIm = sys.C(1,1)*exc.harmonic.r + ...
      2*(1-sys.epsilon_a)*exc.harmonic.r*Omega*1i*h;

% Complex residual
complex_residual_Re = sRe.*Phi_hat +...

complex_residual_Im = sIm.*Phi_hat + ...

% Complex residual
complex_stiffness = sys.C(1,1)*r+2*1i*(r*Omega*(1-sys.epsilon_a)*h-...
    sys.kappa_c*sin(theta_k_mod*h)*sin(theta_k_0));

% Complex residual
complex_residual = complex_stiffness.*Phi_hat +...
        8*sys.epsilon_a*r^2/pi^2 * F_hat_imag - eye(H+1,1);

% Split into real and imaginary component
residual = [real(complex_residual(1:(H+1))); ...
            imag(complex_residual(2:(H+1)))];

end