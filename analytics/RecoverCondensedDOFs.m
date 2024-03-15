function [Q,theta,expDelta] = RecoverCondensedDOFs(sys,exc,r,xi,disorder)
%RECOVERCONDENSEDDOFS Compute the Fourier Coefficients of remaining
% linear DOFs that are ignored in Condensation Equation
%   r - [1, N] vector of frequencies
%   xi - [1, N] vector of clearance normalized amplitudes
%   angle - [1, N] phase angle of nonlinear sector

r  = permute(r,[1 3 2]);
xi  = permute(xi,[1 3 2]);

% Auxilliary Variable
rho = (2/pi) * (1-sys.eN) / (1+sys.eN);

% Triangle Wave Amplitudes
theta = (1+sqrt((1+rho^2)*xi.^2 - rho^2))/(1+rho^2);

% Phase lag of absorber - exp(-1i * Delta)
expDelta = (theta-1)./xi - 1i*rho*theta./xi;

% Fourier Coefficient of 1:1 resonance
Pi = zeros(sys.N_s,1,length(r));
Pi(1,1,:) = -8*sys.epsilon_a*....
    theta.*expDelta.*r.^2 / pi^2;

% Linear transfer function matrix
switch disorder
    case 'tuned'
        % Transfer function matrix with removed absorbers
        % Start in modal basis
        H = pageinv(-r.^2 .* sys.mu + 1i*r.*sys.beta + sys.kappa);
        % Linear displacement with removed absorbers
        Q_lin = pagemtimes(H,exc.F_mod);
        Q_lin = pagemtimes(sys.Phi,Q_lin);
        % Transform back to phsysical coordinates
        H = pagemtimes(sys.Phi,H);
        H = pagemtimes(H,transpose(sys.Phi));
        % Phase of nonlinear DOF
        Phase = angle(Q_lin(1,1,:)./(xi + H(1,1,:).*Pi(1,1,:)));
        % Fourier Coefficient contact forces
        Pi(1,1,:) = sys.Gamma(1)*Pi(1,1,:).*exp(1i*Phase);
    case 'mistuned'
        % Transfer function matrix with removed absorbers
        H = pageinv(-r.^2 .* sys.M_mt + 1i*r.*sys.C + sys.K_mt);
        % Linear displacement with removed absorbers
        Q_lin = pagemtimes(H,exc.F);
        % Phase of nonlinear DOF
        Phase = angle(Q_lin(1,1,:)./(xi + H(1,1,:).*Pi(1,1,:)));
        % Fourier Coefficient contact forces
        Pi(1,1,:) = sys.Gamma_mt(1)*Pi(1,1,:).*exp(1i*Phase);
    otherwise
        error('Case not defined.')
end

% Solve full HB system
Q = squeeze(Q_lin - pagemtimes(H,Pi));

end

