function [F_hat_real,F_hat_imag]=...
    EvaluateSlowForceCoefficient(Phi_hat,N,H,sys)

% Evaluate the real an imaginary part of the slow nonlinear force
% in the time domain without scaling of 8*epsilon_a*r^2 / pi^2
%
% Phi_hat - Fourier coefficients of amplitude envelope
% N - Sampling points per period for AFT
% H - Harmonic Order


% Slow flow in time domain
Phi = FrequencyTime(Phi_hat,N,'Freq_to_Time');

% Auxilliary Variable
rho = (2/pi) * (1-sys.eN) / (1+sys.eN);

% Check if slow flow falls of SIM
i_not_on_SIM = Phi<(sys.Gamma(1)* rho/sqrt(1+rho^2));

% If not: Correct to avoid error but give warning
if any(i_not_on_SIM)
    warning('Solution not on SIM.')
    Phi(i_not_on_SIM) = 1.0001*sys.Gamma(1)* rho/sqrt(1+rho^2);
end

% Triangle Wave Amplitude of absorber
qahat = (sys.Gamma(1)+sqrt((1+rho^2)*Phi.^2 ...
    - (sys.Gamma(1)*rho)^2))/(1+rho^2);

% Real and imaginary part of slow force coefficient
force_real = qahat.*(qahat-sys.Gamma(1))./Phi;
force_imag = rho*(qahat.^2)./Phi;

% Fourier coefficients
F_hat_real = FrequencyTime(force_real,H,'Time_to_Freq');
F_hat_imag = FrequencyTime(force_imag,H,'Time_to_Freq');

end