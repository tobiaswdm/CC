function [force_real,force_imag]=...
    EvaluateSlowForceCoefficient(slow_flow,sys)

% Evaluate the real an imaginary part of the slow nonlinear force
% in the time domain without scaling of 8*epsilon_a*r^2 / pi^2
%
% slow_flow - real amplitude envelope in time

% Auxilliary Variable
rho = (2/pi) * (1-sys.eN) / (1+sys.eN);

% Check if slow flow falls of SIM
i_not_on_SIM = slow_flow<(sys.Gamma(1)* rho/sqrt(1+rho^2));

% Correct to avoid error but give warning
if any(i_not_on_SIM)
    warning('Solution not on SIM.')
    slow_flow(i_not_on_SIM) = 1.0001*sys.Gamma(1)* rho/sqrt(1+rho^2);
end

% Triangle Wave Amplitude of absorber
qahat = (sys.Gamma(1)+sqrt((1+rho^2)*slow_flow.^2 ...
    - (sys.Gamma(1)*rho)^2))/(1+rho^2);


force_real = qahat.*(qahat-sys.Gamma(1))./slow_flow;
force_imag = rho*(qahat.^2)./slow_flow;

end