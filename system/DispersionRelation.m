function [r,r_removed] = DispersionRelation(k,sys)
    
r = sqrt(1+4*sys.kappa_c*sin(k*pi/sys.N_s).^2);
r_removed = r/sqrt(1-sys.epsilon_a);
            
end