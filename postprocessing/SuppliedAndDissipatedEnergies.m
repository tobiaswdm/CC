function [Esupp,Ediss_host,Ediss_VI] = SuppliedAndDissipatedEnergies( ...
    U,UA,TAU,disorder,exc,sys,sol)

% SuppliedAndDissipatedEnergies Compute the total supplied (Esupp) and 
% dissipated energies by the linear dmaping of the host strucutre 
% (Ediss_host) and the instantaneous damping of the VIs (Ediss_VI)

% U - Physical Velocities Oscillators [N_s,N_tau]
% UA - Physical Velocities VI-NESs [N_s, N_tau]
% TAU - Time vector [1, N_tau]

% Supplied energies
Esupp = sum(U.*real(exc.F*exp(1i*exc.harmonic.r*TAU))*sol.dtau,"all");

% Initialize VI
Ediss_VI = 0;

% Dissipated
switch disorder
    case 'tuned'
        % by host
        Ediss_host = sum(U.*(sys.C*U),'all')*sol.dtau;

        % by VI-NESs
        for i = 2:length(TAU)
            % Detect sectors with impacts
            impact = UA(:,i-1) ~= UA(:,i);
            
            if any(~isnan(impact))
                Ediss_VI = Ediss_VI + ...
                           0.5*U(impact,i-1)'*sys.M(impact,impact)*U(impact,i-1)+...
                           0.5*UA(impact,i-1)'*sys.Ma(impact,impact)*UA(impact,i-1)-...
                           0.5*U(impact,i)'*sys.M(impact,impact)*U(impact,i)-...
                           0.5*UA(impact,i)'*sys.Ma(impact,impact)*UA(impact,i);
            end
        end

    case 'mistuned'
        % by host
        Ediss_host = sum(U.*(sys.C_mt*U),'all')*sol.dtau;

        % by VI-NESs
        for i = 2:length(TAU)
            % Detect sectors with impacts
            impact = UA(:,i-1) ~= UA(:,i);
            
            if any(~isnan(impact))
                Ediss_VI = Ediss_VI + ...
                           0.5*U(impact,i-1)'*sys.M_mt(impact,impact)*U(impact,i-1)+...
                           0.5*UA(impact,i-1)'*sys.Ma(impact,impact)*UA(impact,i-1)-...
                           0.5*U(impact,i)'*sys.M_mt(impact,impact)*U(impact,i)-...
                           0.5*UA(impact,i)'*sys.Ma(impact,impact)*UA(impact,i);
            end
        end
    otherwise
        error('Case not defined.')
end


end

