function [sys,exc] = BuildSystem(sys,exc,disorder)
    
    %% Tuned system
    switch disorder
        case 'tuned' % Build tuned system
            
            % Modal wavenumbers
            sys.k = 0:floor(sys.N_s/2);

            % Eigenfrequencies
            [sys.r_k,sys.r_k_noabs,~] = DispersionRelation(sys.k,sys);

            % Build tuned stiffness matrix
            if sys.N_s>2
                K = (1+2*sys.kappa_c)*eye(sys.N_s);
                K = K - diag(sys.kappa_c*ones(1,sys.N_s-1),1) ...
                    - diag(sys.kappa_c*ones(1,sys.N_s-1),-1);
                K(1,sys.N_s) =  - sys.kappa_c;
                K(sys.N_s,1) =  - sys.kappa_c;
                sys.K = K;
            elseif sys.N_s == 2
                K = (1+sys.kappa_c)*eye(sys.N_s);
                K(1,sys.N_s) =  - sys.kappa_c;
                K(sys.N_s,1) =  - sys.kappa_c;
            elseif sys.N_s == 1
                K = 1;
            end

        
            % Build tuned mass matrices

            % Removed abosrbers (use in nonlinear sim)
            sys.M = (1-sys.epsilon_a)*eye(sys.N_s);
            % Fixed absorbers
            sys.M_fixed = eye(sys.N_s);
            sys.mu = (1-sys.epsilon_a)*eye(sys.N_s); % Modal damping matrix removed abosrbers (use in nonlinear sim)
            sys.mu_fixed = eye(sys.N_s); % Modal damping matrix fixed absorbers
            
            
            % Determine linear modes shapes and eigenfrequencies    
            % Differ between even and uneven number of sectors N_s
            sys.k_sort = zeros(1,sys.N_s); % Match system matrices to wavenumber
            if rem(sys.N_s,2)
                k = 1;
                j = (0:(sys.N_s-1))';
                r_eigen = ones(1,sys.N_s);
                Phi = ones(sys.N_s);
                Phi(:,1) = Phi(:,1)/sqrt(sys.N_s);
                sys.k_sort(1) = 0;
                    for i = 2:2:sys.N_s
                        sys.k_sort(i) = k;
                        sys.k_sort(i+1) = k;
                        r_eigen(i) = sys.r_k(k+1);
                        r_eigen(i+1) = sys.r_k(k+1);
                        Phi(:,i) = cos(2*pi*j*k/sys.N_s);
                        Phi(:,i+1) = sin(2*pi*j*k/sys.N_s);
                        Phi(:,i) = Phi(:,i)/norm(Phi(:,i));
                        Phi(:,i+1) = Phi(:,i+1)/norm(Phi(:,i+1));
                        k=k+1;
                    end
            else
                k = 1;
                j = (0:(sys.N_s-1))';
                r_eigen = ones(1,sys.N_s);
                Phi = ones(sys.N_s);
                Phi(:,1) = Phi(:,1)/sqrt(sys.N_s);
                sys.k_sort(1) = 0;
                for i = 2:2:(sys.N_s-1)
                    sys.k_sort(i) = k;
                    sys.k_sort(i+1) = k;
                    r_eigen(i) = sys.r_k(k+1);
                    r_eigen(i+1) = sys.r_k(k+1);
                    Phi(:,i) = cos(2*pi*j*k/sys.N_s);
                    Phi(:,i+1) = sin(2*pi*j*k/sys.N_s);
                    Phi(:,i) = Phi(:,i)/norm(Phi(:,i));
                    Phi(:,i+1) = Phi(:,i+1)/norm(Phi(:,i+1));
                    k=k+1;
                end
                r_eigen(sys.N_s) = sys.r_k(end);
                Phi(:,sys.N_s) = (-1).^j / sqrt(sys.N_s);
                sys.k_sort(sys.N_s) = sys.N_s/2;
            end   
            sys.Phi = Phi; % Modal transformation matrix
            

            % Modal stiffness matrix
            sys.kappa = diag(r_eigen.^2);

            % Damping matrices
            sys.beta = diag(2*sys.D*r_eigen);   % Modal 
            sys.C = Phi*sys.beta*transpose(Phi); % Physical 
   

            % Absorber mass matrix
            sys.Ma = sys.epsilon_a * eye(sys.N_s);
            
            % Reference resonance amplitude
            sys.qref = (2*sys.D*sys.r_k(exc.k+1)^2)^-1;
            
            % Clearance
            sys.Gamma = sys.Gamma_Scale*sys.qref*ones(2*sys.N_s,1);

            % Force directions
            sys.WN = zeros(sys.N_s,2*sys.N_s);
            sys.WN_mod = zeros(sys.N_s,2*sys.N_s);
            sys.WNA = zeros(sys.N_s,2*sys.N_s);
            for i = 1:sys.N_s
                % Modal
                sys.WN_mod(:,(2*i-1):(2*i)) = [-Phi(i,:);Phi(i,:)]';
                % Absorber
                sys.WNA(i,(2*i-1):(2*i)) = [1 -1];
                % Physical
                sys.WN(i,(2*i-1):(2*i)) = [-1 1];
            end

            % Complex excitation vector
            exc.F = exp(2*pi*1i*exc.k/sys.N_s*(0:(sys.N_s-1))');
            exc.F_mod = transpose(Phi)*exc.F;

            sys.tuned_param_set = true;

        case 'mistuned'
            if ~isfield(sys,'tuned_param_set')
                error('Compute tuned system matrices first!')
            end

            % Local eigenfrequency mistuning
            if sys.sigma_omega ~= 0 

                % Take samples and correct to zero mean
                sys.delta_omega = sys.sigma_omega * randn(sys.N_s,1);
                sys.delta_omega = sys.delta_omega - mean(sys.delta_omega);

                % Convert to grounding stiffness mistuning
                sys.delta_kg = (1+sys.delta_omega).^2 - 1;
                
                % Mistuned matrix
                sys.K_mt = sys.K + diag(sys.delta_kg);
                sys.kappa_mt = transpose(sys.Phi)*sys.K_mt*sys.Phi;
            else
                sys.K_mt = sys.K;
                sys.kappa_mt = sys.kappa;
                sys.delta_omega = zeros(sys.N_s,1);
                sys.delta_kg = zeros(sys.N_s,1);
            end

            % Local clearance mistuning
            if sys.sigma_g~= 0 

                % Take samples and correct to zero mean
                sys.delta_g = sys.sigma_g * randn(sys.N_s,1);
                sys.delta_g = sys.delta_g - mean(sys.delta_g);     
                
                % Update clearances
                sys.Gamma_mt = zeros(2*sys.N_s,1);
                sys.Gamma_mt(1:2:(2*sys.N_s)) = ...
                    sys.Gamma(1:2:(2*sys.N_s)).*(1+sys.delta_g);
                sys.Gamma_mt(2:2:(2*sys.N_s)) = ...
                    sys.Gamma(2:2:(2*sys.N_s)).*(1+sys.delta_g);
            else
                sys.Gamma_mt = sys.Gamma;
                sys.delta_g = zeros(sys.N_s,1);
            end

            % Mass matrices
            sys.mu_mt = sys.mu;
            sys.M_mt = sys.M;
            sys.M_fixed_mt = sys.M_fixed;
            
            % Remove or fix absorber in first sector by adjusting
            % Mass matrix and increase clearance in first sector
            if isfield(sys,'absorber_malfunction')
                switch sys.absorber_malfunction
                    case 'fixed'
                        sys.M_mt(1,1) = sys.M_mt(1,1)+sys.epsilon_a;
                        sys.mu_mt = transpose(sys.Phi)*sys.M_mt*sys.Phi;
                        sys.Gamma_mt(1) = 1e6 * sys.qref;
                        sys.Gamma_mt(2) = 1e6 * sys.qref;
                    case 'removed'
                        sys.Gamma_mt(1) = 1e6 * sys.qref;
                        sys.Gamma_mt(2) = 1e6 * sys.qref;
                        sys.M_fixed_mt(1) = sys.M_fixed_mt(1)-...
                            sys.epsilon_a;
                    otherwise
                        error('Case not defined.')
                end

            end
            

            % Compute new eigenfrequencies
            % Fixed absorbers
            [Phi_mt,r_k_mt_sq] = eig(sys.K_mt,sys.M_fixed_mt);
            [sys.r_k_mt,isort] = ...
                sort(transpose(sqrt(diag(r_k_mt_sq))));
     
            % Removed absorbers
            sys.r_k_noabs_mt = ...
                sort(transpose(sqrt(eig(sys.K_mt,sys.M_mt))));

            % Damping matrix
            if sys.adjustC
                % Mass-normalized Eigenvectors of system with fixed absorbers
                Phi_mt = Phi_mt(:,isort);
                Phi_mt = Phi_mt./repmat(sqrt( ...
                    diag(Phi_mt'*sys.M_fixed_mt*Phi_mt)'),[sys.N_s,1]);
                
                % Inverse
                Phi_mt_inv = Phi_mt\eye(sys.N_s);               
                
                % Damping Matrix in Physical Coordinates
                sys.C_mt = Phi_mt_inv'*diag(2*sys.r_k_mt*sys.D)* Phi_mt_inv;
                % Damping Matrix in Modal Coordinates of Tuned system
                sys.beta_mt = sys.Phi' * sys.C_mt * sys.Phi;
            else
                sys.beta_mt = sys.beta;
                sys.C_mt = sys.C;
            end

        case 'mistuned_defined'
            if ~isfield(sys,'tuned_param_set')
                error('Compute tuned system matrices first!')
            end

            % Local eigenfrequency mistuning
            if isfield(sys,'delta_omega') 

                % Convert to grounding stiffness mistuning
                sys.delta_kg = (1+sys.delta_omega).^2 - 1;
                
                % Mistuned matrix
                sys.K_mt = sys.K + diag(sys.delta_kg);
                sys.kappa_mt = transpose(sys.Phi)*sys.K_mt*sys.Phi;

            else
                sys.K_mt = sys.K;
                sys.kappa_mt = sys.kappa;
                sys.delta_omega = zeros(sys.N_s,1);
                sys.delta_kg = zeros(sys.N_s,1);
            end

            % Local clearance mistuning
            if isfield(sys,'delta_g')     
                
                % Update clearances
                sys.Gamma_mt = zeros(2*sys.N_s,1);
                sys.Gamma_mt(1:2:(2*sys.N_s)) = ...
                    sys.Gamma(1:2:(2*sys.N_s)).*(1+sys.delta_g);
                sys.Gamma_mt(2:2:(2*sys.N_s)) = ...
                    sys.Gamma(2:2:(2*sys.N_s)).*(1+sys.delta_g);

            else

                sys.Gamma_mt = sys.Gamma;
                sys.delta_g = zeros(sys.N_s,1);
                
            end

            % Mass matrices
            sys.mu_mt = sys.mu;
            sys.M_mt = sys.M;
            sys.M_fixed_mt = sys.M_fixed;
            
            % Remove or fix absorber in first sector by adjusting
            % Mass matrix and increase clearance in first sector
            if isfield(sys,'absorber_malfunction')
                switch sys.absorber_malfunction
                    case 'fixed'
                        sys.M_mt(1,1) = sys.M_mt(1,1)+sys.epsilon_a;
                        sys.mu_mt = transpose(sys.Phi)*sys.M_mt*sys.Phi;
                        sys.Gamma_mt(1) = 1e6 * sys.qref;
                        sys.Gamma_mt(2) = 1e6 * sys.qref;
                    case 'removed'
                        sys.Gamma_mt(1) = 1e6 * sys.qref;
                        sys.Gamma_mt(2) = 1e6 * sys.qref;
                        sys.M_fixed_mt(1) = sys.M_fixed_mt(1)-...
                            sys.epsilon_a;
                    otherwise
                        error('Case not defined.')
                end

            end
            
            % Compute new eigenfrequencies
            % Fixed absorbers
            [Phi_mt,r_k_mt_sq] = eig(sys.K_mt,sys.M_fixed_mt);
            [sys.r_k_mt,isort] = ...
                sort(transpose(sqrt(diag(r_k_mt_sq))));
     
            % Removed absorbers
            sys.r_k_noabs_mt = ...
                sort(transpose(sqrt(eig(sys.K_mt,sys.M_mt))));

            % Damping matrix
            if sys.adjustC
                % Mass-normalized Eigenvectors of system with fixed absorbers
                Phi_mt = Phi_mt(:,isort);
                Phi_mt = Phi_mt./repmat(sqrt( ...
                    diag(Phi_mt'*sys.M_fixed_mt*Phi_mt)'),[sys.N_s,1]);
                
                % Inverse
                Phi_mt_inv = Phi_mt\eye(sys.N_s);               
                
                % Damping Matrix in Physical Coordinates
                sys.C_mt = Phi_mt_inv'*diag(2*sys.r_k_mt*sys.D)* Phi_mt_inv;
                % Damping Matrix in Modal Coordinates of Tuned system
                sys.beta_mt = sys.Phi' * sys.C_mt * sys.Phi;
            else
                sys.beta_mt = sys.beta;
                sys.C_mt = sys.C;
            end

        otherwise
            error('Case not defined')

    end
   
    
end