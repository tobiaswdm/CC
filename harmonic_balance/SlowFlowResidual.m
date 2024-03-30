function [residual] = SlowFlowResidual(Phi_hat_ReIm,H,k_mod,sys,exc)
%SLOWFLOWRESIDUAL Get the Harmonic Balance residual for the 
% almost quasi periodic response (AQPR)
%
% Phi_hat - 

% Modulation wavenumber phase angle
theta_k_mod =2*pi*k_mod/sys.N_s;

% Excitation wavenumber phase angle
theta_k_0 =2*pi*exc.k/sys.N_s;

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
[force_real,force_imag]=...
    EvaluateSlowForceCoefficient(slow_flow,sys);

% Get Fourier coefficients of dynamic components force coefficients
% without scaling of 8*epsilon_a*r^2 / pi^2

% Fundamental of real component
F_hat_real_0 = mean(force_real);

% All harmonics of imaginary component
F_hat_imag = FrequencyTime(force_imag,H,'Time_to_Freq');


% Evaluate nonlinear dispersion relation and group velocity to get
% fast excitaiton frequency r and the slow modulation frequency Omega
[r,Omega] = ...
    NonlinearDispersionRelation(Phi_hat_0,F_hat_real_0,sys,exc);

% Build residual equation

% Harmonics
h = (0:H)';

% Complex residual
complex_stiffness = sys.C(1,1)+2*1i*(r*Omega*(1-sys.epsilon_a)*h-...
    sys.kappa_c*sin(theta_k_mod*h)*sin(theta_k_0));

% Complex residual
complex_residual = complex_stiffness.*Phi_hat +...
        8*sys.epsilon_a*r^2/pi^2 * F_hat_imag - eye(H+1,1);

% Split into real and imaginary component
residual = [real(complex_residual(1:(H+1))); ...
            imag(complex_residual(2:(H+1)))];

end

function [force_real,force_imag]=...
    EvaluateSlowForceCoefficient(slow_flow,sys)

% Evaluate the real an imaginary part of the slow nonlinear force
% in the time domain without scaling of 8*epsilon_a*r^2 / pi^2

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
force_imag = rho*qahat./slow_flow;

end

function [r,Omega] = ...
    NonlinearDispersionRelation(Phi_hat_0,F_hat_real_0,sys,exc)
% Evaluate nonlinear dispersion relation to get r and Omega

    
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

