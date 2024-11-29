function [qhat_practically_stable,qhat_stable,qhat_unstable,r, ...
    qhat_unstable_synchloss, qhat_unstable_modulation] = ...
    StabilityAnalysisGSR(c,sys,sol,exc)
%STABILITYANALYSIS Study stability along contour plot
%
% c - low level contour estimation

% Find beginning of level curves in c
c(2,floor(c(2,:))==c(2,:)) = NaN;

% Amplitudes in reference sector
xi = c(2,2:end);

% Frequencies
r = c(1,2:end);

% Initialize amplitudes
qhat_practically_stable = nan(1,length(r));
qhat_stable = nan(1,length(r));
qhat_unstable = nan(1,length(r));

% Criteria
qhat_unstable_synchloss = nan(1,length(r));
qhat_unstable_modulation = nan(1,length(r));

parfor (i = 1:length(r), sol.N_Workers)
    
    sol_loop = sol;
    exc_loop = exc;
    sys_loop = sys;

    % Initialize excitaiton structure
    exc_loop.harmonic.r = r(i);
    
    if (sys_loop.kappa_c == 0) || ...
       (exc_loop.k == 0 || exc_loop.k == sys_loop.N_s/2)

        % Choose half of the group velocity period or at least 300 periods
        sol_loop.N_Tau = 300;

    else

        % Evaluate as many periods as determined by group velocity
        [~,~,R_g] = DispersionRelation(exc_loop.k,sys_loop);

        % Choose half of the group velocity period or at least 300 periods
        sol_loop.N_Tau = max(ceil(exc_loop.harmonic.r/R_g/2),300);

    end

    
    if ~isnan(xi(i))
        
        %
        [stable,qhat] = SlowFlowStability('GSR', ...
                sys_loop, ...
                exc_loop, ...
                xi(i),r(i),'tuned');

        if stable
            qhat_stable(i) = qhat(1);

            qhat_ana = xi(i)*sys.Gamma(1);
            % Get initial conditions
            [sol_loop.q0,sol_loop.u0,sol_loop.qa0,sol_loop.ua0] = ...
            GSRInititalConditions(qhat_ana,sys_loop,exc_loop);
    
            % Configure integrator
            sol_loop = ConfigureIntegrator(sol_loop,sys_loop,exc_loop,...
                'no change',false,'tuned');

            % Transient of 1000 Periods
            sol_loop.NP_trans = sol_loop.N_P*1000; 
    
            % Simulation
            [ETA,~,~,UA,~] = ...
                MoreauIntegration(sys_loop,exc_loop,sol_loop,'tuned');
    
            % Extract Amplitudes and Significant impacts
            [qhat,qhatstd] = ...
                MeanAmplitude(sys_loop.Phi*ETA(:,2:end),sol_loop.N_Sample);
            N_SIPP = CountImpacts(UA,sol_loop.N_Tau,'average');
    
            % Relative deviation from analytical model
            xi_dev = abs((mean(qhat)/sys_loop.Gamma(1)-xi(i))/xi(i)); 
       
    
            % Check if solution diverged
            if all(N_SIPP >= 1.99) && ... % All sectors in 1:1 resonance
                all(qhatstd./qhat <= 5e-2) % Amplitude Modulating?
                
                qhat_practically_stable(i) = mean(qhat);
                
                qhat_stable(i) = NaN;
    
            else
    
                % Return which criteria were responsible
                % for stability loss
                if any(N_SIPP < 1.99)
                    qhat_unstable_synchloss(i) = qhat_ana;
                end
    
                if any(qhatstd./qhat > 5e-2)
                    qhat_unstable_modulation(i) = qhat_ana;
                end
    
            end

        else
            qhat_unstable(i) = qhat(1);
        end

    end
             
end

end










