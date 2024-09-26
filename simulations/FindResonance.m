function [Q_mean, Q_std, N_SIPP, varargout] = ...
                     FindResonance(sys,sol,exc,r_steps,disorder)

L = length(r_steps);
Q_mean = zeros(1,L);
Q_std = zeros(1,L);
N_SIPP = zeros(sys.N_s,L);

if strcmp(disorder,'mistuned')
    LF = zeros(1,L);
    IPR = zeros(1,L);
end

% Set first excitation frequency
exc.harmonic.r = r_steps(1);
sol = ConfigureIntegrator(sol,sys,exc,'random',true,disorder);
sol.N_tau = 0;
sol.N_Save = 1;

% Simulate through transient and extract IC
[sol.q0,sol.qa0,sol.u0,sol.ua0,~] = MoreauIntegration(sys,exc,sol,disorder);

for j = 1:L
    exc.harmonic.r = r_steps(j);
    
    sol = ConfigureIntegrator(sol,sys,exc,'no change',false,disorder);
    
    % Additional waittime for absorber malfunction
    % due to light damping of linear respnse in single sector
    if isfield(sys,'absorber_malfunction')
        sol.NP_trans = 300*sol.N_P;
    end
    
    % Simulation
    [ETA,QA,CHI,UA,~] = MoreauIntegration(sys,exc,sol,disorder);
    
    % Extract performance measures
    N_SIPP(:,j) = CountImpacts(UA,sol.N_Tau,'average');
    [ampmean,ampstd] = MeanAmplitude(sys.Phi*ETA(:,2:end),sol.N_Sample);
    switch disorder
        case 'tuned'
            Q_mean(j) = mean(ampmean);
            Q_std(j) = sqrt(mean(ampstd.^2));
        case 'mistuned'
            [Q_mean(j),kmax] = max(ampmean);
            Q_std(j) = ampstd(kmax);
            [IPR(j),LF(j)] = LocalizationMeasures(ampmean,sys);
        otherwise
            error('Case not defined.')
    end
    
    % Assign new starting conditions
    sol.q0 = ETA(:,end);
    sol.u0 = CHI(:,end);
    sol.qa0 = QA(:,end);
    sol.ua0 = UA(:,end);
end

[Q_mean,i] = max(Q_mean,[],2);
Q_std = Q_std(i);

switch disorder
    case 'tuned'
        N_SIPP = mean(N_SIPP(:,i));
    case 'mistuned'
        N_SIPP = N_SIPP(:,i);
        varargout{1} = IPR(i);
        varargout{2} = LF(i);
    otherwise
        error('Case not defined.')
end

end