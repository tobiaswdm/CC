function [r,Omega,Phi_hat] = ...
    SolveSlowFlowEquation(Phi_hat_0,k_mod,sys,exc)
%SOLVESLOWFLOWEQUATION Solve Slow Flow Equton with HB
%
% Phi_hat_0 - Initial guess of HB solution
% k_mod - Modulation wavenumber

% Harmonic order
H = length(Phi_hat_0)-1;

% Split into real and imaginary component
Phi_hat_init = [real(Phi_hat_0(1:(H+1))); ...
                imag(Phi_hat_0(2:(H+1)))];

% Solve using fsolve
Phi_hat_ReIm = fsolve(@(x) SlowFlowResidual(x,H,k_mod,sys,exc),...
                      Phi_hat_init);

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
[force_real,~]=...
    EvaluateSlowForceCoefficient(slow_flow,sys);

% Get Fourier coefficients of dynamic components force coefficients
% without scaling of 8*epsilon_a*r^2 / pi^2
% Fundamental of real component
F_hat_real_0 = mean(force_real);

% Get fast and slow frequency
[r,Omega] = ...
    NonlinearDispersionRelation(Phi_hat_ReIm(1),F_hat_real_0,sys,exc);
end

