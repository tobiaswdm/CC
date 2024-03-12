%% Build the system

% Tuned
[sys,exc] = BuildSystem(sys,exc,'tuned');

% Draw Dispersion Diagram
DrawDispersion(sys,color.reference,savepath);

% Mistuned system
[sys,exc] = BuildSystem(sys,exc, ...
            simsetup.TimeSimulation.disorder);

%% Time integration

% Configure Integrator settings
sol = ConfigureIntegrator(sol,sys,exc,...
    simsetup.TimeSimulation.initial_conditions, ...
    simsetup.TimeSimulation.cut_transient,...
    'tuned'); % Tuned system
sol_mt = ConfigureIntegrator(sol,sys,exc,...
    simsetup.TimeSimulation.initial_conditions, ...
    simsetup.TimeSimulation.cut_transient,...
    'mistuned'); % Mistuned system

% Simulation mistuned system
[ETA_mt,QA_mt,~,~,~] = MoreauIntegration(sys,exc,sol_mt,'mistuned');
Q_mt = sys.Phi * ETA_mt;

% Simulation tuned system
[ETA,QA,~,~,TAU] = MoreauIntegration(sys,exc,sol,'tuned');
Q = sys.Phi * ETA;

%% Post-processing

% Amplitudes - mistuned system
[qhat_mt,qhat_mt_std] = MeanAmplitude(Q_mt(:,2:end),sol.N_Sample);

% Amplitudes - tuned system
[qhat,qhat_std] = MeanAmplitude(Q(:,2:end),sol.N_Sample);