function [F_hat]=EvaluateSlowForceCoefficient(Phi_hat,N,H,sys)

% Evaluate the slow nonlinear force
% in the time domain without scaling of -8*epsilon_a*r^2 / pi^2
%
% Phi_hat - Fourier coefficients of amplitude envelope
% N - Sampling points per period for AFT
% H - Harmonic Order


% Slow flow in time domain
Phi = FrequencyTime(Phi_hat,N,'Freq_to_Time');

% Auxilliary Variable
rho = (2/pi) * (1-sys.eN) / (1+sys.eN);

% Check if slow flow falls of SIM
i_not_on_SIM = abs(Phi)<(sys.Gamma(1)* rho/sqrt(1+rho^2));

% If not: Correct to avoid error but give warning
if any(i_not_on_SIM)
    warning('Solution not on SIM.')
    Phi(i_not_on_SIM) = 1.0001*sys.Gamma(1)* rho/sqrt(1+rho^2)*...
                        exp(1i*angle(Phi(i_not_on_SIM)));
end

% Triangle Wave Amplitude of absorber
qahat = (sys.Gamma(1)+sqrt((1+rho^2)* abs(Phi).^2 ...
        - (sys.Gamma(1)*rho)^2))/(1+rho^2);

% Slow Force Coefficient
force = qahat.*(qahat*(1-1i*rho)-sys.Gamma(1))./abs(Phi) .* ...
        exp(1i*angle(Phi));

% Fourier coefficients
F_hat = FrequencyTime(force,H,'Time_to_Freq');

end