function [Phi_hat,Omega] = ...
    SolveSlowFlowEquation(Phi_hat_0,Omega_0,k_mod,sys,exc)
%SOLVESLOWFLOWEQUATION Solve Slow Flow Equton with HB
%
% Phi_hat_0 - Initial guess of HB solution
% k_mod - Modulation wavenumber

% Harmonic order
H = (length(Phi_hat_0)-1)/2;

% Split into real and imaginary component
Phi_hat_init = [real(Phi_hat_0); imag(Phi_hat_0);1000*Omega_0];
% Delete real part of first harmonic
Phi_hat_init(H+2) = [];

% Solve using fsolve
options = optimoptions('fsolve','Display','iter');
Phi_hat_ReIm = fsolve(@(x) SlowFlowResidual(x,H,k_mod,sys,exc),...
                      Phi_hat_init,options);


Omega = Phi_hat_ReIm(end)/1000;
% Convert to complex Fourier Coefficients
Phi_hat = [Phi_hat_ReIm(1:(H+1));10;Phi_hat_ReIm((H+2):(2*H))]...
+1i*Phi_hat_ReIm((2*H)+(1:(2*H+1)));

% Excitation wavenumber phase angle
%theta_k_0 = 2*pi*exc.k/sys.N_s;

% Evaluate nonlinear dispersion relation for average amplitude level
% get slow modulation frequency Omega (group velocity)
%[~,Omega] = NonlinearDispersionRelation(abs(Phi_hat(H+1)),theta_k_0,sys);

end

