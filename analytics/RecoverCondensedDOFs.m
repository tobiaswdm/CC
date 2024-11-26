function [Q,theta,expDelta] = RecoverCondensedDOFs(sys,exc,r,xi,pattern, ...
    varargin)
%RECOVERCONDENSEDDOFS Compute the Fourier Coefficients of remaining
% linear DOFs that are ignored in Condensation Equation
%   r - [1, N] vector of frequencies
%   xi - [1, N] vector of clearance normalized amplitudes
%   angle - [1, N] phase angle of nonlinear sector
%   pattern - 'single' or 'opposing'
%   varagin{1} contains disorder 'tuned' or 'mistuned'

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

% Linear transfer function matrix
switch pattern
    case 'single'
        Pi(1,1,:) = -8*sys.epsilon_a*....
        theta.*expDelta.*r.^2 / pi^2;

        switch varargin{1}
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
                H = pageinv(-r.^2 .* sys.M_mt + 1i*r.*sys.C_mt + sys.K_mt);
                % Linear displacement with removed absorbers
                Q_lin = pagemtimes(H,exc.F);
                % Phase of nonlinear DOF
                Phase = angle(Q_lin(1,1,:)./(xi + H(1,1,:).*Pi(1,1,:)));
                % Fourier Coefficient contact forces
                Pi(1,1,:) = sys.Gamma_mt(1)*Pi(1,1,:).*exp(1i*Phase);
            otherwise
                error('Case not defined.')
        end
    case 'opposing'
        
        Pi(1,1,:) = -8*sys.epsilon_a*....
                    theta.*expDelta.*r.^2 / pi^2;

        % Build vector of relevant wavenumbers
        k = zeros(1,sys.N_s);
        if rem(sys.N_s,2)~=0
            error('FRS not defined for even number of sectors.')
        else
            k(2:2:sys.N_s) = sys.k(2:end);
            k(3:2:(sys.N_s-1)) = sys.k(2:end-1);
        end

        % Consider only even or odd modal wavenumbers
        if rem(exc.k,2)~=0
            k = rem(k,2) ~= 0;    
        else
            k = rem(k,2) == 0;  
        end

        % Linear transfer function matrix
        H = pageinv(-r.^2 .* sys.mu(k,k) + 1i*r.*sys.beta(k,k) ...
            + sys.kappa(k,k));

        % Linear displacement with removed absorbers
        Q_lin = pagemtimes(H,exc.F_mod(k));
        Q_lin = pagemtimes(sys.Phi(:,k),Q_lin);

        % Transform back to phsysical coordinates
        H = pagemtimes(sys.Phi(:,k),H);
        H = pagemtimes(H,transpose(sys.Phi(:,k)));

        % Phase of nonlinear DOFs
        Phase = angle(Q_lin(1,1,:)./(xi + 2*H(1,1,:).*Pi(1,1,:)));
        
        % Fourier Coefficient contact forces
        Pi(1,1,:) = sys.Gamma(1)*Pi(1,1,:).*exp(1i*Phase);

        if rem(exc.k,2)~=0
            Pi(sys.N_s/2+1,1,:) = -Pi(1,1,:);
        else
            Pi(sys.N_s/2+1,1,:) = Pi(1,1,:);
        end

    otherwise
        error('Case not defined.')
end

% Solve full HB system
Q = squeeze(Q_lin - pagemtimes(H,Pi));

end

