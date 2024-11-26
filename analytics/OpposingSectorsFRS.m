function [Gamma_Scale,Xi,R] = OpposingSectorsFRS(xi,r,sys,exc)
% Determine the frequency amplitude surface (FRS) of the synchronization in
% two opposing sectors for the tuned configuration
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

% Build vector of wavenumbers
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
H = pageinv(-r_perm.^2 .* sys.mu(k,k) + 1i*r_perm.*sys.beta(k,k) ...
    + sys.kappa(k,k));
Hnn = pagemtimes(sys.Phi(1,k),H);
Hnn = repmat(squeeze(pagemtimes(Hnn,transpose(sys.Phi(1,k)))),[1,Nxi]);
Q_nn_lin = repmat(transpose((-(1-sys.epsilon_a)*r.^2 +...
        2*sys.D*1i*sys.r_k(exc.k+1)*r + ....
        sys.r_k(exc.k+1)^2).^-1),[1,Nxi]);

% FRS
Gamma_Scale = abs(Q_nn_lin./(Xi + 2*Hnn.*Pi))/sys.qref; % Clearance


end

