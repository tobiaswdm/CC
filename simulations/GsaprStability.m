%% Build tuned system

% Tuned
[sys,exc] = BuildSystem(sys,exc,'tuned');

% Draw Dispersion Diagram
DrawDispersion(sys,color,savepath);

%% Compute ESIM

% Auxilliary Variable
rho = (2/pi) * (1-sys.eN) / (1+sys.eN);
% Minimum amplitude - turning point of SIM
% Including safety of 0.1% higher amplitude
xi_min = 1.001 * rho/sqrt(1+rho^2);

% Clearance normalized amplitudes
xi = logspace(log10(xi_min),...
    log10(simsetup.GsaprStability.xi_max),...
    simsetup.GsaprStability.Nxi);
r = linspace(simsetup.GsaprStability.r_range(1),...
    simsetup.GsaprStability.r_range(2), ...
    simsetup.GsaprStability.Nr);

% Linear FRF
q_fixed = abs(ComputeLinearResponse(r,sys,exc,'tuned','fixed_absorbers'));
q_fixed = q_fixed(1,:);
q_removed = abs(ComputeLinearResponse(r,sys,exc,'tuned','removed_absorbers'));
q_removed = q_removed(1,:);

% ESIM
[Gamma_Scale,Xi,R] = AllSectorsESIM(xi,r,sys,exc);

% Backbone
theta = (1+sqrt((1+rho^2)*xi.^2 - rho^2))/(1+rho^2);
varpi_bb = sqrt(xi./((1-sys.epsilon_a)*xi + ...
    8*sys.epsilon_a*theta.*((theta-1)./xi)/pi^2));
Delta = acos((theta-1)./xi);
Gamma_scale_bb = 2*sys.D./abs((-(1-sys.epsilon_a)*varpi_bb.^2 + ...
    2*sys.D*1i*varpi_bb + 1).*xi - ...
    8*sys.epsilon_a*theta.*varpi_bb.^2 .* exp(-1i*Delta)/pi^2);


% Plot ESIM
figure(2);
surf(R/sys.r_k(exc.k+1),Xi,Gamma_Scale,'EdgeColor','none',...
    'DisplayName','ESIM')
hold on;
plot3(varpi_bb,xi,Gamma_scale_bb,'--k','LineWidth',3,...
    'DisplayName','Backbone')
title('ESIM GSAPR')
box on;
legend;
colormap turbo
xlabel('$\varpi$')
ylabel('$\xi$')
zlabel('$\Gamma/\hat{q}_\mathrm{ref}$')
set(gca,'YScale','log')
axis tight;

% Plot Level Curves
figure(3);
contour(R/sys.r_k(exc.k+1),Xi,Gamma_Scale,10,'LineWidth',1.5,...
    'DisplayName','ESIM')
hold on;
plot(varpi_bb,xi,'--k','LineWidth',3,'DisplayName','Backbone')
legend;
h=colorbar;
h.Label.Interpreter = 'latex';
h.Label.String = "$\Gamma/\hat{q}_\mathrm{ref}$";
title('ESIM GSAPR')
box on;
colormap turbo
xlabel('$\varpi$')
ylabel('$\xi$')
set(gca,'YScale','log')
axis tight;

% Get Level curves at clearance
c = contourc(r,xi,Gamma_Scale',sys.Gamma_Scale*[1 1]);

% Determine max amplitude FRF
r_plot = c(1,2:end);
qhat = c(2,2:end)*sys.Gamma(1);
qhat(floor(c(2,:))==c(2,:)) = NaN;

% Coarsen contour for stability analysis
c = CoarsenContour(c,...
    simsetup.GsaprStability.stepsize);

% Study asymptotic and practical stability of tuned system
[qhat_stable,qhat_unstable,r_num,qhat_unstable_synchloss, ...
qhat_unstable_amplitudedev,qhat_unstable_modulation] = ...
StabilityAnalysisGsapr(c,sys,sol,exc);

% Stability plot
figure(4);
hold on;
plot(r/sys.r_k(exc.k+1),q_fixed/sys.qref,...
    'LineWidth',.5,'Color',color.reference,'DisplayName', ...
    'Fixed abs.')
plot(r/sys.r_k(exc.k+1),q_removed/sys.qref,'-.',...
    'LineWidth',.5,'Color',color.reference,'DisplayName', ...
    'Removed abs.')
plot(r_plot/sys.r_k(exc.k+1),qhat/sys.qref,...
            'LineWidth',1.5,'Color',color.ies,'DisplayName', ...
            'GSAPR')
scatter(r_num/sys.r_k(exc.k+1),qhat_unstable/sys.qref,20,'MarkerFaceColor',color.show,...
    'MarkerEdgeColor','k','Displayname','Unstable')
scatter(r_num/sys.r_k(exc.k+1),qhat_stable/sys.qref,20,'MarkerFaceColor',...
    myColors('green'),'MarkerEdgeColor','k','Displayname','Stable')
set(gca,'YScale','log')
axis tight;
box on;
xlabel('$\varpi$')
ylabel('$\hat{q}/\hat{q}_\mathrm{ref}$')
legend;
savefig([savepath 'frequency_amplitude_stability.fig'])

% Check criterion
figure(5)
tiledlayout(1,3);
nexttile;
hold on;
plot(r/sys.r_k(exc.k+1),q_fixed/sys.qref,...
    'LineWidth',.5,'Color',color.reference,'HandleVisibility', ...
    'off')
plot(r/sys.r_k(exc.k+1),q_removed/sys.qref,'-.',...
    'LineWidth',.5,'Color',color.reference,'HandleVisibility', ...
    'off')
plot(r_plot/sys.r_k(exc.k+1),qhat/sys.qref,...
            'LineWidth',1.5,'Color',color.ies,'HandleVisibility', ...
    'off')
scatter(r_num(isnan(qhat_unstable_synchloss))/sys.r_k(exc.k+1), ...
    qhat_unstable(isnan(qhat_unstable_synchloss))/sys.qref,20, ...
    'MarkerFaceColor',color.show,'MarkerEdgeColor','k','HandleVisibility',...
    'off', 'MarkerEdgeAlpha',.25,'MarkerFaceAlpha',.25)
scatter(r_num/sys.r_k(exc.k+1),qhat_stable/sys.qref,20,'MarkerFaceColor',...
    myColors('green'),'MarkerEdgeColor','k','HandleVisibility','off', ...
    'MarkerEdgeAlpha',.25,'MarkerFaceAlpha',.25)
scatter(r_num/sys.r_k(exc.k+1), qhat_unstable_synchloss/sys.qref,20, ...
    'MarkerFaceColor',myColors('magenta'),'MarkerEdgeColor','k', ...
    'Displayname', 'Criterion Broken')
set(gca,'YScale','log')
axis tight;
box on;
xlabel('$\varpi$')
ylabel('$\hat{q}/\hat{q}_\mathrm{ref}$')
title('Criterion 1.')
legend;
nexttile;
hold on;
plot(r/sys.r_k(exc.k+1),q_fixed/sys.qref,...
    'LineWidth',.5,'Color',color.reference,'HandleVisibility', ...
    'off')
plot(r/sys.r_k(exc.k+1),q_removed/sys.qref,'-.',...
    'LineWidth',.5,'Color',color.reference,'HandleVisibility', ...
    'off')
plot(r_plot/sys.r_k(exc.k+1),qhat/sys.qref,...
            'LineWidth',1.5,'Color',color.ies,'HandleVisibility', ...
    'off')
scatter(r_num(isnan(qhat_unstable_amplitudedev))/sys.r_k(exc.k+1), ...
    qhat_unstable(isnan(qhat_unstable_amplitudedev))/sys.qref,20, ...
    'MarkerFaceColor',color.show,'MarkerEdgeColor','k','Displayname', ...
    'Unstable', 'MarkerEdgeAlpha',.25,'MarkerFaceAlpha',.25)
scatter(r_num/sys.r_k(exc.k+1),qhat_stable/sys.qref,20,'MarkerFaceColor',...
    myColors('green'),'MarkerEdgeColor','k','Displayname','Stable', ...
    'MarkerEdgeAlpha',.25,'MarkerFaceAlpha',.25)
scatter(r_num/sys.r_k(exc.k+1), qhat_unstable_amplitudedev/sys.qref,20, ...
    'MarkerFaceColor',myColors('magenta'),'MarkerEdgeColor','k', ...
    'Displayname', 'Criterion Broken')
set(gca,'YScale','log')
axis tight;
box on;
xlabel('$\varpi$')
ylabel('$\hat{q}/\hat{q}_\mathrm{ref}$')
title('Criterion 2.')
nexttile;
hold on;
plot(r/sys.r_k(exc.k+1),q_fixed/sys.qref,...
    'LineWidth',.5,'Color',color.reference,'DisplayName', ...
    'Fixed abs.')
plot(r/sys.r_k(exc.k+1),q_removed/sys.qref,'-.',...
    'LineWidth',.5,'Color',color.reference,'DisplayName', ...
    'Removed abs.')
plot(r_plot/sys.r_k(exc.k+1),qhat/sys.qref,...
            'LineWidth',1.5,'Color',color.ies,'DisplayName', ...
            'GSAPR')
scatter(r_num(isnan(qhat_unstable_modulation))/sys.r_k(exc.k+1), ...
    qhat_unstable(isnan(qhat_unstable_modulation))/sys.qref,20, ...
    'MarkerFaceColor',color.show,'MarkerEdgeColor','k','Displayname', ...
    'Unstable', 'MarkerEdgeAlpha',.25,'MarkerFaceAlpha',.25)
scatter(r_num/sys.r_k(exc.k+1),qhat_stable/sys.qref,20,'MarkerFaceColor',...
    myColors('green'),'MarkerEdgeColor','k','Displayname','Stable', ...
    'MarkerEdgeAlpha',.25,'MarkerFaceAlpha',.25)
scatter(r_num/sys.r_k(exc.k+1), qhat_unstable_modulation/sys.qref,20, ...
    'MarkerFaceColor',myColors('magenta'),'MarkerEdgeColor','k', ...
    'Displayname', 'Criterion Broken')
set(gca,'YScale','log')
axis tight;
box on;
xlabel('$\varpi$')
ylabel('$\hat{q}/\hat{q}_\mathrm{ref}$')
title('Criterion 3.')
savefig([savepath 'stability_criteria.fig'])
