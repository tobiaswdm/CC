function [F_hat_real,F_hat_imag]=...
    EvaluateSlowCouplingForceCoefficient(Phi_hat,Gamma_hat,theta_k_mod,...
                                        theta_k_0,N,H)

% Evaluate the real an imaginary part of the slow coupling force
% in the time domain without scaling of kappa_c
%
% Phi_hat - Fourier coefficients of amplitude envelope
% Gamma_hat - Fourier coefficients of phase envelope
% theta_k_mod - Modulation inter sector phase angle
% theta_k_0 - Excitation inter-sector phase angle
% N - Sampling points per period for AFT
% H - Harmonic Order

% Coefficients of neighboring sectors
h = (0:H)'; % Harmonics
Phi_hat_minus = Phi_hat .* exp(-1i*h*theta_k_mod);
Phi_hat_plus = Phi_hat .* exp(1i*h*theta_k_mod);
Gamma_hat_minus = Gamma_hat .* exp(-1i*h*theta_k_mod);
Gamma_hat_plus = Gamma_hat .* exp(1i*h*theta_k_mod);

% Transform to time domain
Phi_minus = FrequencyTime(Phi_hat_minus,N,'Freq_to_Time');
Phi_plus = FrequencyTime(Phi_hat_plus,N,'Freq_to_Time');
Gamma_minus = FrequencyTime(Gamma_hat_minus,N,'Freq_to_Time');
Gamma_plus = FrequencyTime(Gamma_hat_plus,N,'Freq_to_Time');
Gamma = FrequencyTime(Gamma_hat,N,'Freq_to_Time');

% Compute force in time domain
force_real = Phi_plus.*cos(Gamma_plus-Gamma+theta_k_0) + ...
             Phi_minus.*cos(Gamma_minus-Gamma-theta_k_0);
force_imag = Phi_plus.*sin(Gamma_plus-Gamma+theta_k_0) + ...
             Phi_minus.*sin(Gamma_minus-Gamma-theta_k_0);

% Fourier coefficients
F_hat_real = FrequencyTime(force_real,H,'Time_to_Freq');
F_hat_imag = FrequencyTime(force_imag,H,'Time_to_Freq');

end

