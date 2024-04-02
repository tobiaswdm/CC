function [r,Omega] = ...
    NonlinearDispersionRelation(Phi_hat_0,theta_k_0,sys)
% Evaluate nonlinear dispersion relation to get fast frequency r 
% and group velocity Omega
%
% Phi_hat_0 - real valued harmonic amplitude
% theta_k0 - Excitation inter-sector phase angle 2*pi*exc.k/sys.N_s

% Auxilliary Variable
rho = (2/pi) * (1-sys.eN) / (1+sys.eN);

% Triangle Wave Amplitudes
qahat = (sys.Gamma(1)+sqrt((1+rho^2)*Phi_hat_0.^2 - ...
        (sys.Gamma(1)*rho)^2))/(1+rho^2);

% Scaling constant of dispersion
c = 1-sys.epsilon_a + (8*sys.epsilon_a/pi^2)*...
    (qahat.*(qahat-sys.Gamma(1))./(Phi_hat_0.^2));

    
% Nonlinear dispersion relation
r = sqrt((1+4*sys.kappa_c*sin(theta_k_0/2)^2)./c);

% Nonlinear group velocity
Omega = 2*pi/sys.N_s * sys.kappa_c * sin(theta_k_0) ./ ...
        sqrt(c .* (1+4*sys.kappa_c*sin(theta_k_0/2)^2));

end

