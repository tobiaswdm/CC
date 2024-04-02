function [F_hat_real,F_hat_imag] =...
    EvaluateSlowExcitationForceCoefficient(Gamma_hat,N,H)
% Evaluate the real an imaginary part of the slow excitation force
% in the time domain without
%
% Gamma_hat - Fourier coefficients of phase envelope
% N - Sampling points per period for AFT
% H - Harmonic Order

Gamma = FrequencyTime(Gamma_hat,N,'Freq_to_Time');

% Compute force in time domain
force_real = -cos(Gamma);
force_imag = sin(Gamma);

% Fourier coefficients
F_hat_real = FrequencyTime(force_real,H,'Time_to_Freq');
F_hat_imag = FrequencyTime(force_imag,H,'Time_to_Freq');

end

