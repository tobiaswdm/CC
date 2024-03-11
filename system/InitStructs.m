%% System parameters
% Tuned system
sys.N_s = 10;
sys.kappa_c = 0;
sys.epsilon_a = 0.02;
sys.D = 1e-3;
sys.Gamma_Scale = 0;    % Gamma = qref*Gamma_scale
sys.eN = 0.8;

% Mistuned system
sys.sigma_omega = 0;
sys.sigma_g = 0;

%% Excitation
exc.type = 'harmonic';
exc.k = 0;
exc.wavedirection = 'forward';
% Harmonic excitation
exc.harmonic.r = 0;
% Sweep
exc.sweep.r0 = 0;   % Start frequency
exc.sweep.re = 0;   % End frequency
exc.sweep.tau = 0;  % Sweep duration
exc.sweep.r_gref = 0; % r to calculate reference amplitude


%% Solver
sol.alpha = 1;          % Only change if necassary
sol.mode = 'smoreau';   % Only change if necassary
sol.solver = 'JOR';     % Only change if necassary
sol.maxiter = 1000;     % Only change if necassary
sol.tol = 1e-4;         % Only change if necassary
sol.N_Tau = 300;        % Number of stationary excitation periods to simulate
sol.N_P = 1000;         % Sample pseudo period with 1000 steps
sol.N_Sample = 1000;    % Sample pseudo period with 1000 steps in post processing
sol.N_Workers = 25;     % Maximum number of workers

%% Simulation Setups

% =========================================================================
% Variation Coupling and Clearance with MCS
% =========================================================================

% Number of Monte Carlos Simulations per nominal configuration
simsetup.VariationCouplingAndClearanceMCS.N_MCS = 50;

% Parameters of nominal coupling
% Range
simsetup.VariationCouplingAndClearanceMCS.Range_kappa_c = [1e-4,3];
% Scaling of discrete points in interval
simsetup.VariationCouplingAndClearanceMCS.Scaling_kappa_c = 'logarithmic';
% Number of discrete points in interval
simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c = 10;

% Parameters of nominal clearance normalized by tuned resonance amplitude
% Range
simsetup.VariationCouplingAndClearanceMCS.Range_GammaScale = [1e-4,3];
% Scaling of discrete points in interval
simsetup.VariationCouplingAndClearanceMCS.Scaling_GammaScale = 'logarithmic';
% Number of discrete points in interval
simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale = 10;

% Frequency stepping range
simsetup.VariationCouplingAndClearanceMCS.r_scale = [0.98 1.03];
% Steps in frequency stepping interval
simsetup.VariationCouplingAndClearanceMCS.N_rSteps = 50;

% =========================================================================
% Variation Coupling and Clearance
% =========================================================================

% Parameters of nominal coupling
% Range
simsetup.VariationCouplingAndClearance.Range_kappa_c = [1e-4,3];
% Scaling of discrete points in interval
simsetup.VariationCouplingAndClearance.Scaling_kappa_c = 'logarithmic';
% Number of discrete points in interval
simsetup.VariationCouplingAndClearance.Number_kappa_c = 10;

% Parameters of nominal clearance normalized by tuned resonance amplitude
% Range
simsetup.VariationCouplingAndClearance.Range_GammaScale = [1e-4,3];
% Scaling of discrete points in interval
simsetup.VariationCouplingAndClearance.Scaling_GammaScale = 'logarithmic';
% Number of discrete points in interval
simsetup.VariationCouplingAndClearance.Number_GammaScale = 10;

% Frequency stepping range
simsetup.VariationCouplingAndClearance.r_scale = [0.98 1.03];
% Steps in frequency stepping interval
simsetup.VariationCouplingAndClearance.N_rSteps = 50;

% =========================================================================
% Convergence of Monte Carlo Simulation
% =========================================================================

% Number of Monte Carlos Simulations per nominal configuration
% Maximum Number of MCSs
simsetup.ConvergenceMCS.N_MCS = 10000;
% Number of MCS Permutations for converegence analysis
simsetup.ConvergenceMCS.Shuffle_MCS = 1000;

% Frequency stepping range
simsetup.ConvergenceMCS.r_scale = [0.98 1.03];
% Steps in frequency stepping interval
simsetup.ConvergenceMCS.N_rSteps = 50;









