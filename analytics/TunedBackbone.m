function [g_scale,q_scale,g_opt] = TunedBackbone(sys,stability)

D = sys.D;
eN = sys.eN;
epsa = sys.epsilon_a;

%% Backbone

rho = (2/pi)*(1-eN)/(1+eN);

xi = logspace(log10((1+1e-5)*rho/sqrt(1+rho^2)),log10(100),1000);
theta = (1+sqrt((1+rho^2)*xi.^2 - rho^2))/(1+rho^2);
Delta = acos((theta-1)./xi);

varpi = sqrt(xi./((1-epsa)*xi + 8*epsa*theta.*cos(Delta)/pi^2));

g_scale = 2*D./abs((-(1-epsa)*varpi.^2 + 2*D*1i*varpi + 1).*xi - ...
                    8*epsa*theta.*varpi.^2 .* exp(-1i*Delta)/pi^2);

[g_opt,ii] = max(g_scale);

switch stability
    case 'stable'
        g_scale = g_scale(ii:end);
        q_scale = xi(ii:end).*g_scale;
    case 'full'
        q_scale = xi.*g_scale;
    otherwise
        error('Case not defined.')
end

end