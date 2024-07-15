function [Gamma_Scale,Xi,R] = SingleSectorFRS(xi,r,sys,exc,disorder)
% Determine the frequency amplitude surface (FRS) of the synchronization in a single sector
% for either the tuned or the mistuned configuration
%   
%   Complex_Phase of the localized sector
%   Gamma_Scale size of the clearance of the localized sector scaled by
%   tuned reference amplitude
%
%   r = [1,Nr] - frequencies to eveluate the FRS at
%   xi = [1,Nxi] - clearance normalized ampltiude to evaluate the FRS at

% Get lengths
Nr = length(r);
Nxi = length(xi);

% Compute grid to evaluate FRS over
[Xi,R] = meshgrid(xi,r);

% Auxilliary Variable
rho = (2/pi) * (1-sys.eN) / (1+sys.eN);

% Triangle Wave Amplitudes
Theta = repmat((1+sqrt((1+rho^2)*xi.^2 - rho^2))/(1+rho^2),[Nr,1]);

% Phase lag of absorber - exp(-1i * Delta)
expDelta = (Theta-1)./Xi - 1i*rho*Theta./Xi;

% Fourier Coefficient of 1:1 resonance
Pi = -8*sys.epsilon_a*Theta.*expDelta.*R.^2 / pi^2;

% Frequencies
r_perm = permute(r,[1 3 2]);

% Linear transfer function matrix
switch disorder
    case 'tuned'
        H = pageinv(-r_perm.^2 .* sys.mu + 1i*r_perm.*sys.beta + sys.kappa);
        Hnn = pagemtimes(sys.Phi(1,:),H);
        Hnn = repmat(squeeze(pagemtimes(Hnn,transpose(sys.Phi(1,:)))),[1,Nxi]);
        Q_nn_lin = repmat(transpose((-(1-sys.epsilon_a)*r.^2 +...
                2*sys.D*1i*sys.r_k(exc.k+1)*r + ....
                sys.r_k(exc.k+1)^2).^-1),[1,Nxi]);
    case 'mistuned'
        H = pageinv(-r_perm.^2 .* sys.M_mt + 1i*r_perm.*sys.C + sys.K_mt);
        Hnn = repmat(squeeze(H(1,1,:)),[1 Nxi]);
        Q_lin = pagemtimes(H(1,:,:),exc.F);
        Q_nn_lin = repmat(squeeze(Q_lin),[1,Nxi]);
    otherwise
        error('Case not defined.')
end

Gamma_Scale = abs(Q_nn_lin./(Xi + Hnn.*Pi))/sys.qref; % Clearance


end

