function [Q] = RecoverCondensedDOFs(sys,exc,r,xi,phase)
%RECOVERCONDENSEDDOFS Compute the Fourier Coefficients of remaining
% linear DOFs that are ignored in Condensation Equation
%   r - [1, N] vector of frequencies
%   xi - [1, N] vector of clearance normalized amplitudes
%   angle - [1, N] phase angle of nonlinear sector

r  = permute(r,[1 3 2]);
xi  = permute(xi,[1 3 2]);
phase = permute(phase,[1 3 2]);

% Auxilliary Variable
rho = (2/pi) * (1-sys.eN) / (1+sys.eN);

% Triangle Wave Amplitudes
theta = (1+sqrt((1+rho^2)*xi.^2 - rho^2))/(1+rho^2);

% Phase lag of absorber - exp(-1i * Delta)
expDelta = cos((theta-1)./xi) - 1i*sin(rho*theta./xi);

% Fourier Coefficient of 1:1 resonance
Pi = zeros(sys.N_s,1,length(r));
Pi(1,1,:) = -8*sys.epsilon_a*sys.qref*sys.Gamma_Scale*exp(1i*phase)....
    .*theta.*expDelta.*r.^2 / pi^2;

% Linear transfer function matrix
switch disorder
    case 'tuned'
        % Start in modal basis
        H = pageinv(-r.^2 .* sys.mu + 1i*r.*sys.beta + sys.kappa);
        Q_lin = pagemtimes(H,exc.F_mod);
        Q_lin = pagemtimes(sys.Phi,Q_lin);
        H = pagemtimes(sys.Phi,H);
        H = pagemtimes(H,transpose(sys.Phi));
    case 'mistuned'
        H = pageinv(-r.^2 .* sys.M_mt + 1i*r.*sys.C + sys.K_mt);
        Q_lin = pagemtimes(H,exc.F);
    otherwise
        error('Case not defined.')
end

Q = squeeze(Q_lin - pagemtimes(H,Pi));

end

