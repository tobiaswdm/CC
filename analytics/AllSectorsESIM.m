function [Gamma_Scale,Xi,R] = AllSectorsESIM(xi,r,sys,exc)
% Determine the extendend SIM (ESIM) of the localization in a single sector
% for either the tuned or the mistuned configuration
%   
%   Complex_Phase of the localized sector
%   Gamma_Scale size of the clearance of the localized sector scaled by
%   tuned reference amplitude
%
%   r = [1,Nr] - frequencies to eveluate the ESIM at
%   xi = [1,Nxi] - clearance normalized ampltiude to evaluate the ESIM at

% Get lengths
Nr = length(r);
Nxi = length(xi);

% Compute grid to evaluate ESIM over
[Xi,R] = meshgrid(xi,r);

% Auxilliary Variable
rho = (2/pi) * (1-sys.eN) / (1+sys.eN);

% Triangle Wave Amplitudes
Theta = repmat((1+sqrt((1+rho^2)*xi.^2 - rho^2))/(1+rho^2),[Nr,1]);

% Phase lag of absorber - exp(-1i * Delta)
expDelta = (Theta-1)./Xi - 1i*rho*Theta./Xi;

% Fourier Coefficient of 1:1 resonance
Pi = -8*sys.epsilon_a*Theta.*expDelta.*R.^2 / pi^2;

% Linear eigenfrequency with fixed absorbers
r_k0 = sys.r_k(exc.k+1);

%Evaluate ESIM
Gamma_Scale = 2*sys.D*r_k0^2 ./ ...
               ((-(1-sys.epsilon_a)*R.^2 + 2*sys.D*1i*R*r_k0 + ...
               r_k0^2)*Xi+Pi);

end

