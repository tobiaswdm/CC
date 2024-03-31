function [r,Omega] = ...
    NonlinearDispersionRelation(Phi_hat_0,F_hat_real_0,sys,exc)
% Evaluate nonlinear dispersion relation to get fast frequency r 
% and group velocity Omega
%
% Phi_hat_0 - real valued harmonic amplitude
% F_hat_real_0 - harmonic amplitude of real nonlinear force component

    
% Excitation inter-sector phase angle
theta_k0 = 2*pi*exc.k/sys.N_s;

% Nonlinear dispersion relation
r = sqrt(Phi_hat_0*(1+4*sys.kappa_c*sin(theta_k0/2)^2)/...
    ((1-sys.epsilon_a)*Phi_hat_0+8*sys.epsilon_a*F_hat_real_0/pi^2));

% Nonlinear group velocity
Omega = 2*pi/sys.N_s * sys.kappa_c * sin(theta_k0) * Phi_hat_0 / ...
    sqrt(Phi_hat_0*(1+4*sys.kappa_c*sin(theta_k0/2)^2) *...
    ((1-sys.epsilon_a)*Phi_hat_0 + 8*sys.epsilon_a*F_hat_real_0/pi^2));

end