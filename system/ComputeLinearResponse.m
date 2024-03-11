function [Q] = ComputeLinearResponse(r,sys,exc,disorder,config)

% Initialize complex amplitudes
Q = zeros(sys.N_s,length(r));

switch config % Fixed or removed absorbers
    case 'fixed_absorbers'

        switch disorder % Tuned or mistuned
            case 'tuned'
                for i = 1:length(r)
                    Q(:,i) = (-r(i)^2*sys.M_fixed + 1i*r(i)*sys.C + sys.K)\exc.F;
                end
            case 'mistuned'
                for i = 1:length(r)
                    Q(:,i) = (-r(i)^2*sys.M_fixed_mt + 1i*r(i)*sys.C + sys.K_mt)\exc.F;
                end
            otherwise
                error('Case not defined.')
        end

    case 'removed_absorbers'

        switch disorder % Tuned or mistuned
            case 'tuned'
                for i = 1:length(r)
                    Q(:,i) = (-r(i)^2*sys.M + 1i*r(i)*sys.C + sys.K)\exc.F;
                end
            case 'mistuned'
                for i = 1:length(r)
                    Q(:,i) = (-r(i)^2*sys.M_mt + 1i*r(i)*sys.C + sys.K_mt)\exc.F;
                end
            otherwise
                error('Case not defined.')
        end

    otherwise
        error('Case not defined.')
end

end