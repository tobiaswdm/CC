function [qhat_stable,qhat_unstable,r] = ...
StabilityAnalysisGsapr(c,sys,sol,exc)
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
qhat_stable = nan(1,length(r));
qhat_unstable = nan(1,length(r));

% Only evaluate final 50 Periods
sol.N_Tau = 50;

parfor (i = 1:length(r), sol.N_Workers)
%for i = 1:length(r)
    
    sol_loop = sol;
    exc_loop = exc;
    sys_loop = sys;

    % Initialize excitaiton structure
    exc_loop.harmonic.r = r(i);
    
    if ~isnan(xi(i))
        qhat_ana = xi(i)*sys.Gamma(1);
        % Get initial conditions
        [sol_loop.q0,sol_loop.u0,sol_loop.qa0,sol_loop.ua0] = ...
        GsaprInititalConditions(qhat_ana,sys_loop,exc_loop);

        % Configure integrator
        sol_loop = ConfigureIntegrator(sol_loop,sys_loop,exc_loop,...
            'no change',false,'tuned');
        % Transient of 1000 Periods
        sol_loop.NP_trans = sol_loop.N_P*1000; 

        % Simulation
        [ETA,~,~,UA,~] = ...
            MoreauIntegration(sys_loop,exc_loop,sol_loop,'tuned');

        % Extract Amplitudes and Significant impacts
        [qhat,~] = ...
            MeanAmplitude(sys_loop.Phi*ETA(:,2:end),sol_loop.N_Sample);
        N_SIPP = CountImpacts(UA,sol_loop.N_Tau,'average');

        % Relative deviation from analytical model
        xi_dev = abs((mean(qhat)/sys_loop.Gamma(1)-xi(i))/xi(i)); 

        % Check if solution diverged
        if all(N_SIPP >= 1.99) && ... % All sectors in 1:1 resonance
            xi_dev<=0.75 % Relative deviation of amplitude max 50%
            
            qhat_stable(i) = mean(qhat);

        else

            qhat_unstable(i) = qhat_ana;

        end

    end
             
end

end










