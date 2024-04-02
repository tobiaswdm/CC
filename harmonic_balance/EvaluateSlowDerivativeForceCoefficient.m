function [F_hat_real] = ...
    EvaluateSlowDerivativeForceCoefficient(Phi_hat,Gamma_hat,N,H)

% Evaluate the real an imaginary part of the slow coupling force
% in the time domain without scaling of kappa_c
%
% Phi_hat - Fourier coefficients of amplitude envelope
% Gamma_hat - Fourier coefficients of phase envelope
% N - Sampling points per period for AFT
% H - Harmonic Order

% Derivative of Gamma
h = (0:H)'; % Harmonics
dGamma_hat = Gamma_hat .* (1i*h);

dGamma = FrequencyTime(dGamma_hat,N,'Freq_to_Time');
Phi = FrequencyTime(Phi_hat,N,'Freq_to_Time');

% Compute force in time domain
force_real = -dGamma.*Phi;

% Fourier coefficients
F_hat_real = FrequencyTime(force_real,H,'Time_to_Freq');

end

