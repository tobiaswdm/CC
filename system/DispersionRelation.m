function [r,r_removed,R_g] = DispersionRelation(k,sys)

% Dispersion relation fixed absorbers
r = sqrt(1+4*sys.kappa_c*sin(k*pi/sys.N_s).^2);

% Transform to removed absorbers
r_removed = r/sqrt(1-sys.epsilon_a);

% Group Velocity
R_g = 2*pi/sys.N_s * sys.kappa_c * sin(2*k*pi/sys.N_s) ./ ...
        sqrt((1+4*sys.kappa_c*sin(k*pi/sys.N_s).^2));
            
end