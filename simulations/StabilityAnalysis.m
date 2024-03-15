function [qhat_practically_stable,qhat_stable,qhat_unstable,r] =...
    StabilityAnalysis(c,sys,sol,exc,disorder,stability,practical_stability)
%STABILITYANALYSIS Study stability along contour plot
%
% c - low level contour estimation
% stability - study asymptotic stability? - true, false
% practical_stability - study practical stability? - true, false

% Determine maximum amplitude allong contour and its kin. accessible ranges
[qhat_max_ana,~,r] = ...
    LocalizedFrequencyAmplitudeCurve(c,sys,exc,disorder);

% Find beginning of level curves in c
c(2,floor(c(2,:))==c(2,:)) = NaN;

% Amplitudes in reference sector
xi = c(2,2:end);

% Initialize amplitudes
qhat_practically_stable = nan(1,length(r));
qhat_stable = nan(1,length(r));
qhat_unstable = nan(1,length(r));

% Only evaluate final 50 Periods
sol.N_Tau = 50;

parfor (i = 1:length(r), sol.N_Workers)  
    
    % Only consider cases where kinematic constraint is fulfilled
    if ~isnan(qhat_max_ana(i))

        sol_loop = sol;
        exc_loop = exc;
        sys_loop = sys;

        % Initialize excitaiton structure
        exc_loop.harmonic.r = r(i);

        if stability
            % Get initial conditions
            [sol_loop.q0,sol_loop.u0,sol_loop.qa0,sol_loop.ua0] = ...
                LocalizedInitialConditions(sys_loop,exc_loop,xi(i),r(i),...
                disorder,'stability');

            % Configure integrator
            sol_loop = ConfigureIntegrator(sol_loop,sys_loop,exc_loop,...
                'no change',false,disorder);
            % Transient of 1000 Periods
            sol_loop.NP_trans = sol_loop.N_P*1000; 

            % Simulation
            [ETA,~,~,UA,~] = ...
                MoreauIntegration(sys_loop,exc_loop,sol_loop,disorder);
            
            % Extract Amplitudes and Significant impacts
            [qhat,~] = ...
                MeanAmplitude(sys_loop.Phi*ETA(:,2:end),sol_loop.N_Sample);
            N_SIPP = CountImpacts(UA,sol_loop.N_Tau);
            
            % Determine relative deviation of amplitude in localized sector
            switch disorder
                case 'tuned'
                    xi_dev = abs((qhat(1)/sys_loop.Gamma(1)- ...
                        xi(i))/xi(i));
                case 'mistuned'
                    xi_dev = abs((qhat(1)/sys_loop.Gamma_mt(1)- ...
                        xi(i))/xi(i));
                otherwise
                    error('Case not defined.')
            end

            % Check if solution diverged
            if N_SIPP(1) >= 1.99 && ... % Sector in 1:1 resonance
                all(N_SIPP(2:end)==0) && ... % No impacts in remaining secs
                xi_dev<=0.5 % Relative deviation of amplitude max 50%
                
                qhat_stable(i) = max(qhat);
                
            end

        end



        if practical_stability
            
            % if stability was studied before, only consider
            % asymptotically stable solutions
            
            if (stability && ~isnan(qhat_stable(i))) || ~stability
                % Get initial conditions
                [sol_loop.q0,sol_loop.u0,sol_loop.qa0,sol_loop.ua0] = ...
                    LocalizedInitialConditions(sys_loop,exc_loop,xi(i),r(i),...
                    disorder,'practical_stability');
    
                % Configure integrator
                sol_loop = ...
                    ConfigureIntegrator(sol_loop,sys_loop,exc_loop,...
                    'no change',false,disorder);
                % Transient of 1000 Periods
                sol_loop.NP_trans = sol_loop.N_P*1000;

                % Simulation
                [ETA,~,~,UA,~] = ...
                    MoreauIntegration(sys_loop,exc_loop,sol_loop,disorder);
                
                % Extract Amplitudes and Significant impacts
                [qhat,~] = ...
                    MeanAmplitude(sys_loop.Phi*ETA(:,2:end),sol_loop.N_Sample);
                N_SIPP = CountImpacts(UA,sol_loop.N_Tau);

                % Determine relative deviation of amplitude in localized sector
                switch disorder
                    case 'tuned'
                        xi_dev = abs((qhat(1)/sys_loop.Gamma(1)- ...
                            xi(i))/xi(i));
                    case 'mistuned'
                        xi_dev = abs((qhat(1)/sys_loop.Gamma_mt(1)- ...
                            xi(i))/xi(i));
                    otherwise
                        error('Case not defined.')
                end

                % Check if solution diverged
                if N_SIPP(1) >= 1.99 && ... % Sector in 1:1 resonance
                    all(N_SIPP(2:end)<1.98) && ... % No 1:1 remaining secs
                    xi_dev<=0.5 % Relative deviation of amplitude max 50%
                    
                    qhat_practically_stable(i) = max(qhat);               
                end

            end

        end
        
        % Assign unstable solution
        if isnan(qhat_practically_stable(i)) && ...
                isnan(qhat_stable(i))
            qhat_unstable(i) = qhat_max_ana(i);
        end

    end

end

end








