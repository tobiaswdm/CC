%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CC (pronounced Sisi) is a tool that performs numerical and/or
% analytical analyses on a Cylic Chain of Oscialltors with Vibro-Impact
% Nonlinear Energy Sinks (VI-NESs)
%
% The Code for CC was written by:
% Tobias Weidemann - (C) 2024
% University of Stuttgart, Germany
% Institute of Aircraft Propulsion Systems
%
% Contact: tobias.weidemann@ila.uni-stuttgart.de
%
% Feel free to use, share and modify under the GPL-3.0 license.
% CC is purely academic and comes with no warranty.
% If you use CC for your own research, please refer to the paper:
%
% T. Weidemann, L. A. Bergman, A. F. Vakakis, M. Krack. (2024)
% "Energy Transfer and Localization in a Forced Cyclic Chain of
% Oscillators with Vibro-Impact Nonlinear Energy Sinks".
% Manuscript submitted to Nonlinear Dynamics
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ETA0,CHI0,QA0,UA0] = GSRInititalConditions(Q,sys,exc)
%LOCALIZEDINITIALCONDITIONS Return initial conditions of HB solution
%
% Q- Complex amplitude vector oscialltors

% Check if only uniform amplitude is given
if length(Q)==1 && sys.N_s ~= 1
    qhat = Q;
    Q = qhat*ones(sys.N_s,1);
end

% Inter-sector phase angle
theta_k0 = 2*pi*exc.k/sys.N_s;

% Auxilliary Variable
rho = (2/pi) * (1-sys.eN) / (1+sys.eN);

% Triangle Wave Amplitude of absorber
qahat = (sys.Gamma(1:sys.N_s)+sqrt((1+rho^2)* abs(Q).^2 ...
    - (sys.Gamma(1:sys.N_s)*rho).^2))/(1+rho^2);

Delta = acos((qahat-sys.Gamma(1:sys.N_s))./abs(Q)); 

% Get modal initital conditions
if length(Q)==1 && sys.N_s ~= 1
    % Get Phase shift if only uniform amplitude is given
    Q = Q./((-(1-sys.epsilon_a)*exc.harmonic.r^2 + 2*sys.D*exc.harmonic.r*...
        sys.r_k(exc.k+1)*1i + sys.r_k(exc.k+1)^2)*abs(Q)-8*sys.epsilon_a*...
        exc.harmonic.r^2*qahat.*exp(-1i*Delta)/pi^2);
end
Q=Q.*exp(1i*theta_k0*(0:(sys.N_s-1))');
% Recall that inv(Phi)=transpose(Phi) for tuned system modes
ETA0 = transpose(sys.Phi)*real(Q);
CHI0 = transpose(sys.Phi)*real(1i*exc.harmonic.r*Q);

% Get initial displacement of absorbers

% Triangle Wave motion
QA0 = 2*qahat.*asin(cos(angle(Q)-Delta))/pi;

% Velocity magnitude
qadot = 2*exc.harmonic.r*qahat/pi;

index = 0<=wrapToPi(angle(Q)-Delta) & ...
    wrapToPi(angle(Q)-Delta)<=pi;

UA0 = zeros(sys.N_s,1); % Velocity
% Assign velocities
UA0(index) = -qadot(index);
UA0(~index) = qadot(~index);


end



