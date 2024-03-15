function [ETA0,CHI0,QA0,UA0] = LocalizedInitialConditions(sys,exc,xi,r,disorder,stability)
%LOCALIZEDINITIALCONDITIONS Return initial conditions of HB solution
%
% xi -scalar clearance normalized amplitude
% r - scalar excitation frequency

[Q,theta,expDelta] = RecoverCondensedDOFs(sys,exc,r,xi,disorder);

% expDelta = exp(-1i*Delta) -> negative sign
Delta = -angle(expDelta);

% Get modal initital conditions
% Recall that inv(Phi)=transpose(Phi) for tuned system modes
ETA0 = transpose(sys.Phi) * Q;
CHI0 = real(1i*r.*ETA0);
ETA0 = real(ETA0);

% Get initial displacement of absorbers

QA0 = zeros(sys.N_s,length(r)); % Displacement
UA0 = zeros(sys.N_s,length(r)); % Velocity


switch disorder
    case 'tuned'
        qahat = sys.Gamma(1)*theta;
    case 'mistuned'
        qahat = sys.Gamma_mt(1)*theta;
    otherwise
        error('Case not defined.')
end

% Triangle Wave motion
QA0(1,:) = 2*qahat*asin(cos(angle(Q(1,:))-Delta))/pi;

% Velocity magnitude
qadot = 2*r.*qahat/pi;

index = 0<=wrapToPi(angle(Q(1,:))-Delta) & ...
    wrapToPi(angle(Q(1,:))-Delta)<=pi;

% Assign velocities
UA0(1,index) = -qadot;
UA0(1,~index) = qadot;

if strcmp(stability,'practical_stability')
    % Assign same initial velocity as oscillator
    % in non localized sector
    UA0(2:end,:) = sys.Phi(2:end,:)*CHI0;

    % Place absorber at cavity wall
    switch disorder
        case 'tuned'
            QA0(2:end) = sys.Gamma(3:2:end) + ...
                sys.Phi(2:end,:)*ETA0;
        case 'mistuned'
            QA0(2:end) = sys.Gamma_mt(3:2:end) + ...
                sys.Phi(2:end,:)*ETA0;
    end
end


end



