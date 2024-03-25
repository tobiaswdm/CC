function [R,K,Xi] = NonlinearDispersion(sys,xi,k)
%NONLINEARDISPERSION Summary of this function goes here
%   Detailed explanation goes here
[K,Xi] = meshgrid(k,xi);

% Auxilliary Variable
rho = (2/pi) * (1-sys.eN) / (1+sys.eN);

% Numerator
Num = 1+2*sys.kappa_c*(1-cos(2*pi*K/sys.N_s));

% Absorber Amplitude
Theta = (1+sqrt((1+rho^2)*Xi.^2 - rho^2))/(1+rho^2);

% cosDelta
cosDelta = (Theta-1)./Xi;

% Denominator
Den = 1-sys.epsilon_a+8*sys.epsilon_a*(Theta./Xi/pi^2).*cosDelta;

R = sqrt(Num./Den);




end

