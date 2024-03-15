%% Build tuned system

% Tuned
[sys,exc] = BuildSystem(sys,exc,'tuned');

% Draw Dispersion Diagram
DrawDispersion(sys,color,savepath);

%% Compute tuned manifold

% Auxilliary Variable
rho = (2/pi) * (1-sys.eN) / (1+sys.eN);
% Minimum amplitude - turning point of SIM
% Including safety of 0.1% higher amplitude
xi_min = 1.001 * rho/sqrt(1+rho^2);

% Clearance normalized amplitudes
xi = logspace(log10(xi_min),...
    log10(simsetup.LocalizationSingleSectorStability.xi_max),...
    simsetup.LocalizationSingleSectorStability.Nxi);
r = linspace(simsetup.LocalizationSingleSectorStability.r_range(1),...
    simsetup.LocalizationSingleSectorStability.r_range(2), ...
    simsetup.LocalizationSingleSectorStability.Nr);

% Plot ESIM
figure(2);
surf(R,Xi,Gamma_Scale,'EdgeAlpha',0)
hold on;
title('Tuned System')
box on;
colormap turbo
xlabel('$r$')
ylabel('$\xi$')
zlabel('$\Gamma/\hat{q}_\mathrm{ref}$')
set(gca,'YScale','log')
axis tight;

% Plot Level Curves
figure(3);
contour(R,Xi,Gamma_Scale,10,'LineWidth',1.5)
hold on;
h=colorbar;
h.Label.Interpreter = 'latex';
h.Label.String = "$\Gamma/\hat{q}_\mathrm{ref}$";
title('Tuned System')
box on;
colormap turbo
xlabel('$r$')
ylabel('$\xi$')
set(gca,'YScale','log')
axis tight;

% Get Level curves at clearance
c = contourc(r,xi,Gamma_Scale,[sys.Gamma_Scale sys.Gamma_Scale]);

% Determine max amplitude FRF
[qhat_max,qhat_max_violated,r_plot] = ...
    LocalizedFrequencyAmplitudeCurve(c,sys,exc,'tuned');

% Coarsen contour for stability analysis
[c_coarse] = CoarsenContour(c,...
    simsetup.LocalizationSingleSectorStability.stepsize);

% Study asymptotic and practical stability of tuned system
[qhat_practically_stable,qhat_stable,qhat_unstable,r] =...
    StabilityAnalysis(c,sys,sol,exc,'tuned',true,true);




