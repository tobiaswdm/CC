%% System parameters
% Tuned system
sys.N_s = 10;               % Number of sectors
sys.kappa_c = 0;            % Linear coupling strength
sys.epsilon_a = 0.02;       % Mass ratio of VI-NES
sys.D = 1e-3;               % Uniform Modal Damping Ratio
sys.Gamma_Scale = 0;        % Clearance normalized by linear resonance amplitude
sys.eN = 0.8;               % Restitution coefficient

% Mistuned system
sys.sigma_omega = 0;        % Rel. Standard deviation of local eigefrequencies
sys.sigma_g = 0;            % Rel. Standard deviation of local clearances
sys.adjustC = false;        % Set true if sys.D should also refer to
                            % modes of the mistuned system and false if
                            % damping matrix of tuned system should also be
                            % used for the mistuned case

%% Excitation
exc.type = 'harmonic';  % Excitation type 'harmonic' or 'sweep'
exc.k = 0;              % Excitation wavenumber
exc.wavedirection = 'forward'; % Direction of traveling wave excitation
% Harmonic excitation
exc.harmonic.r = 1;     % Excitation Frequency
% Sweep
exc.sweep.r0 = 0;   % Start frequency
exc.sweep.re = 0;   % End frequency
exc.sweep.tau = 0;  % Sweep duration


%% Solver
sol.alpha = 1;          % Scaling Parameter of prox solver in Moreau Code
sol.solver = 'JOR';     % Solver type of inclusions ('JOR' or 'fsolve')
sol.maxiter = 1000;     % Maximum number of prox iterations in Moreau sol
sol.tol = 1e-4;         % Tolerance of prox solver
sol.N_Tau = 300;        % Number of stationary excitation periods to simulate
sol.N_P = 1000;         % Sample pseudo period with 1000 steps
sol.N_Sample = 1000;    % Sample pseudo period with 1000 steps in post processing
sol.N_Workers = 1000;   % Maximum number of workers for paralles tasks

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

% =========================================================================
% Single Time Simulation
% =========================================================================

% Type of mistuning - 'tuned', 'mistuned_defined' or 'mistuned'
simsetup.TimeSimulation.disorder = 'mistuned';
% Cut transient response in the beginning?
simsetup.TimeSimulation.cut_transient = true;
% Initial conditions - 'random', 'zero' or 'localized'
simsetup.TimeSimulation.initial_conditions = 'random';

% =========================================================================
% Analytical Study of Synchronization in single sector
% =========================================================================

% Maximum of clearance normalized amplitude
simsetup.SynchronizationSingleSectorAnalytical.xi_max = 20;
% Number of samples of clearance normalized amplitude
simsetup.SynchronizationSingleSectorAnalytical.Nxi = 1000;
% Range of excitaiton frequencies
simsetup.SynchronizationSingleSectorAnalytical.r_range = [0.9 1.1];
% Number of samples of excitation frequencies
simsetup.SynchronizationSingleSectorAnalytical.Nr = 1000;
% Number of Mistuning realizations for frequency amplitude curves
simsetup.SynchronizationSingleSectorAnalytical.N_MCS = 10;

% =========================================================================
% Numerical stability analysis of Synchronization in single sector
% =========================================================================

% Maximum of clearance normalized amplitude
simsetup.SynchronizationSingleSectorStability.xi_max = 20;
% Number of samples of clearance normalized amplitude
simsetup.SynchronizationSingleSectorStability.Nxi = 1000;
% Range of excitaiton frequencies
simsetup.SynchronizationSingleSectorStability.r_range = [0.9 1.1];
% Number of samples of excitation frequencies
simsetup.SynchronizationSingleSectorStability.Nr = 1000;
% Number of Mistuning realizations for frequency amplitude curves
simsetup.SynchronizationSingleSectorStability.N_MCS = 10;
% Take every stepsize-th point of contour for stability analysis
simsetup.SynchronizationSingleSectorStability.stepsize = 200;
% Frequency range for estimation of tuned resonance
simsetup.SynchronizationSingleSectorStability.r_scale = [0.98 1.03];
% Steps in frequency stepping interval
simsetup.SynchronizationSingleSectorStability.N_rSteps = 50;
% Localization Measure 'LF' or 'IPR'
simsetup.SynchronizationSingleSectorStability.LocalizationMeasure = 'LF';

% =========================================================================
% Numerical stability analysis of GSR
% =========================================================================

% Maximum of clearance normalized amplitude
simsetup.GSRStability.xi_max = 20;
% Number of samples of clearance normalized amplitude
simsetup.GSRStability.Nxi = 1000;
% Range of excitaiton frequencies
simsetup.GSRStability.r_range = [0.9 1.1];
% Number of samples of excitation frequencies
simsetup.GSRStability.Nr = 1000;
% Take every stepsize-th point of contour for stability analysis
simsetup.GSRStability.stepsize = 200;

% =========================================================================
% Linear Mistuning Analysis using MCS
% =========================================================================

% Number of Monte Carlos Simulations per nominal configuration
simsetup.LinearMistuningAnalysis.N_MCS = 50;

% Parameters of nominal coupling
% Range
simsetup.LinearMistuningAnalysis.Range_kappa_c = [1e-4,3];
% Scaling of discrete points in interval
simsetup.LinearMistuningAnalysis.Scaling_kappa_c = 'logarithmic';
% Number of discrete points in interval
simsetup.LinearMistuningAnalysis.Number_kappa_c = 10;

% =========================================================================
% Follow the first bifurcation point of the FRS of the LSR
% =========================================================================

% Parameters of clearance-normalized amplitude
% Maximum of clearance normalized amplitude
simsetup.FollowBifurcationFRSofLSR.xi_max = 20;
% Number of samples of clearance normalized amplitude
simsetup.FollowBifurcationFRSofLSR.Nxi = 1000;
% Step Size for stability of frequency amplitude curve
simsetup.FollowBifurcationFRSofLSR.stepsize = 50;

% Parameters of excitation frequency
% Sampling range
simsetup.FollowBifurcationFRSofLSR.r_range = [0.9 1.1];
% Number of samples
simsetup.FollowBifurcationFRSofLSR.Nr = 2000;

% Parameters of nominal coupling
% Range
simsetup.FollowBifurcationFRSofLSR.kappac_range = [1e-4,3];
% Number of discrete points in interval
simsetup.FollowBifurcationFRSofLSR.Nkappac = 10;

% =========================================================================
% Numerical stability analysis of Synchronization in opposing sectors
% =========================================================================

% Maximum of clearance normalized amplitude
simsetup.SynchronizationOpposingSectorsStability.xi_max = 20;
% Number of samples of clearance normalized amplitude
simsetup.SynchronizationOpposingSectorsStability.Nxi = 1000;
% Range of excitaiton frequencies
simsetup.SynchronizationOpposingSectorsStability.r_range = [0.9 1.1];
% Number of samples of excitation frequencies
simsetup.SynchronizationOpposingSectorsStability.Nr = 1000;
% Take every stepsize-th point of contour for stability analysis
simsetup.SynchronizationOpposingSectorsStability.stepsize = 200;

% =========================================================================
% Stability Analysis along Backbone of GSR
% =========================================================================

% Range of inter-sector coupling strengths
simsetup.BackBoneStability.kappac_range = [10^-3, 10^-1];
% Number of inter-sector coupling strengths
simsetup.BackBoneStability.Nkappac = 20;
% Number of clearance-normlaized amplitudes
simsetup.BackBoneStability.Nxi = 30;