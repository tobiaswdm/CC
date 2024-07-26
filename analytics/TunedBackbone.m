function [Gamma_scale,q_scale,Gamma_opt] = TunedBackbone(sys,stability)
% Get Ad-Hoc Analytical backbone curve of the GSR
%
% stability - 'stable' or 'full' ('unstable' not implemented)

D = sys.D;
eN = sys.eN;
epsa = sys.epsilon_a;

%% Backbone

% Auxilliary Variable
rho = (2/pi)*(1-eN)/(1+eN);

% Range of Clearance Normalized Amplitude
xi = logspace(log10((1+1e-5)*rho/sqrt(1+rho^2)),log10(100),1000);

% Clearance Normalized Absorber Amplitude
theta = (1+sqrt((1+rho^2)*xi.^2 - rho^2))/(1+rho^2);
Delta = acos((theta-1)./xi);

% Normalized Frequencies on backbone Curve
varpi = sqrt(xi./((1-epsa)*xi + 8*epsa*theta.*cos(Delta)/pi^2));

% Clearance normalized by linear resonance amplitude on backbone
Gamma_scale = 2*D./abs((-(1-epsa)*varpi.^2 + 2*D*1i*varpi + 1).*xi - ...
                    8*epsa*theta.*varpi.^2 .* exp(-1i*Delta)/pi^2);

% Maximum clearance on backbone (optimum)
[Gamma_opt,ii] = max(Gamma_scale);

% Determine whether you want full backbone our only stable part
switch stability
    case 'stable'
        Gamma_scale = Gamma_scale(ii:end);
        q_scale = xi(ii:end).*Gamma_scale;
    case 'full'
        q_scale = xi.*Gamma_scale;
    otherwise
        error('Case not defined.')
end

end