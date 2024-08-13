function [E] = SpatialEnergies(sys,sol,Q,U,UA,disorder)
%SPATIALENERGIES Compute the Average spatial energies
% over each period of oscillation
%   Q - Displacemet in physical coordinates [N_s,N_tau]
%   U - Velocity in physical coordinates [N_s,N_tau]
%   UA - Velocity of Absorbers in physical coordinates [N_s,N_tau]

movRMS = dsp.MovingRMS(sol.N_Sample);

Q_diff = movRMS(diff([Q;Q(1,:)],[],1)').^2;

switch disorder
    case 'tuned'
        E = 0.5*(...
            movRMS(Q').^2 + ...
            (1-sys.epsilon_a)*movRMS(U').^2 + ...
            0.5*sys.kappa_c*(...
            Q_diff + ...
            circshift(Q_diff,1,2) ...
            ) + ...
            sys.epsilon_a*movRMS(UA').^2 ...
            )';
    case 'mistuned'
        E = 0.5*(...
            (1+sys.delta_kg').*movRMS(Q').^2 + ...
            (1-sys.epsilon_a)*movRMS(U').^2 + ...
            0.5*sys.kappa_c*(...
            Q_diff + ...
            circshift(Q_diff,1,2) ...
            ) + ...
            sys.epsilon_a*movRMS(UA').^2 ...
            )'; 
    otherwise
        error('Case not defined.')
end

end

