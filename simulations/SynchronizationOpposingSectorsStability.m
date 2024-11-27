%% Build tuned system

% Tuned
[sys,exc] = BuildSystem(sys,exc,'tuned');

% Draw Dispersion Diagram
DrawDispersion(sys,color,savepath);

% Number of LSRs
[N_LSR_tuned,N_LSR_mistuned] = NumberOfLSRs(sys.N_s);
disp(['There are ' num2str(N_LSR_tuned) ' coexisiting LSRs in the tuned ' ...
    'system and ' num2str(N_LSR_mistuned) ' coexisiting LSRs in the ' ...
    'mistuned system.'])

%% Compute tuned FRS

% Auxilliary Variable
rho = (2/pi) * (1-sys.eN) / (1+sys.eN);
% Minimum amplitude - turning point of SIM
% Including safety of 0.1% higher amplitude
xi_min = 1.001 * rho/sqrt(1+rho^2);

% Clearance normalized amplitudes
xi = logspace(log10(xi_min),...
    log10(simsetup.SynchronizationOpposingSectorsStability.xi_max),...
    simsetup.SynchronizationOpposingSectorsStability.Nxi);
r = linspace(simsetup.SynchronizationOpposingSectorsStability.r_range(1),...
    simsetup.SynchronizationOpposingSectorsStability.r_range(2), ...
    simsetup.SynchronizationOpposingSectorsStability.Nr);

% Linear FRF
q_fixed = abs(ComputeLinearResponse(r,sys,exc,'tuned','fixed_absorbers'));
q_fixed = q_fixed(1,:);
q_removed = abs(ComputeLinearResponse(r,sys,exc,'tuned','removed_absorbers'));
q_removed = q_removed(1,:);

% FRS
[Gamma_Scale,Xi,R] = OpposingSectorsFRS(xi,r,sys,exc);

% Plot FRS
figure(2);
hold on; box on;
surf(R,Xi,Gamma_Scale,'LineStyle','none')
contour3(R,Xi,Gamma_Scale,4,'LineWidth',1,'Color',color.analytics)
hold off;
box on;
colormap turbo
xlabel('$r$')
ylabel('$\xi$')
zlabel('$\Gamma/\hat{q}_\mathrm{ref}$')
set(gca,'YScale','log')
axis tight;
h=colorbar;
h.Label.Interpreter = 'latex';
h.Label.String = "$\Gamma/\hat{q}_\mathrm{ref}$";
view([0 90])

% Get Level curves at clearance
c = contourc(r,xi,Gamma_Scale',[sys.Gamma_Scale sys.Gamma_Scale]);

% Determine max amplitude FRF
[qhat_max,qhat_max_violated,qhat_syn,r_plot] = ...
    LocalizedFrequencyAmplitudeCurve(c,sys,exc,'opposing','tuned');

% Coarsen contour for stability analysis
c = CoarsenContour(c,...
    simsetup.SynchronizationOpposingSectorsStability.stepsize);

% Study local asymptotic and practical stability 
[qhat_practically_stable,qhat_stable,qhat_unstable, ... 
qhatsynch_practically_stable,qhatsynch_stable,qhatsynch_unstable, ...
IPR,LF,r_num] = StabilityAnalysisLSR(c,sys,sol,exc, ...
                                            'opposing','tuned',true,true);

figure(3);
hold on; box on;
plot(r,q_fixed/sys.qref,...
    'LineWidth',.5,'Color',color.reference,'DisplayName', ...
    'Fixed VI-NESs')
plot(r,q_removed/sys.qref,'-.',...
    'LineWidth',.5,'Color',color.reference,'DisplayName', ...
    'Removed VI-NESs')
plot(r_plot,qhat_syn/sys.qref,'--',...
        'LineWidth',1.5,'Color',color.background,'DisplayName', ...
        'Syn. sector')
plot(r_plot,qhat_max/sys.qref,...
        'LineWidth',1.5,'Color',color.ies,'DisplayName', ...
        'Max. amp.')
plot(r_plot,qhat_max_violated/sys.qref,':',...
            'LineWidth',1.5,'Color',color.show,'DisplayName', ...
            'Viol. kin. constr.')
scatter(r_num,qhatsynch_unstable/sys.qref,20,'MarkerFaceColor',color.show,...
            'MarkerEdgeColor','k','Displayname','Unstable')
scatter(r_num,qhatsynch_stable/sys.qref,20,'MarkerFaceColor',myColors('cyan'),...
            'MarkerEdgeColor','k','Displayname','L. A. Stable')
scatter(r_num,qhatsynch_practically_stable/sys.qref,40,'pentagram',...
            'MarkerFaceColor',myColors('green'),'MarkerEdgeColor','k',...
            'Displayname','Pract. Stable')
hold off;
legend;
set(gca,'YScale','log')
axis tight;
xlabel('$r$')
ylabel('$\hat{q}/\hat{q}_\mathrm{ref}$')
savefig([savepath 'frequency_amplitude_stability.fig'])


if any(~isnan(qhat_practically_stable))
    
    % Extract highest practically stable amplitude
    [xi_max,i_max] = min(qhat_practically_stable/sys.Gamma(1));
    exc.harmonic.r = r_num(i_max);

    % Get initial conditions for time simulation
    % Get initial conditions
    sol.N_Tau = 300;
    [sol.q0,sol.u0,sol.qa0,sol.ua0] = ...
        LocalizedInitialConditions(sys,exc,xi_max,exc.harmonic.r,...
        'opposing','tuned','practical_stability');
    sol = ConfigureIntegrator(sol,sys,exc,...
                    'no change',false,'tuned');
    sol.NP_trans = sol.N_P*1000;

    % Simulation
    [ETA,QA,Chi,~,TAU] = MoreauIntegration(sys,exc,sol,'tuned');
    Q = sys.Phi*ETA;
    U = sys.Phi*Chi;

    % Extract amplitudes
    [qhat,~] = MeanAmplitude(Q(:,2:end),sol.N_Sample);
    

    % Average modal energies tuned system
    [E_mod,E_mod_avg] = ModalEnergies(sys,sol,ETA,Chi);

    figure(4)
    bar(0:floor(sys.N_s/2),E_mod_avg/sum(E_mod_avg),...
        'FaceColor',color.ies)
    hold on;
    xline(exc.k,'--k','LineWidth',1,'Alpha',1)
    xlim([-0.499 floor(sys.N_s/2)+0.499])
    ylim([0 1])
    box on;
    xlabel('Wavenumber - $k$')
    ylabel(['$E_{k,\,\mathrm{avg}}^\mathrm{mod} / ' ...
        'E_{\mathrm{tot, \, avg}}^\mathrm{mod}$'])
    title('Average modal energies - tuned system')
    set(gca,'YScale','log')
    savefig([savepath 'modal_energies_example.fig'])
    
    % Phase space tuned system
    % Get analytical HB approximation
    [Qana,~,~] = ...
        RecoverCondensedDOFs(sys,exc,exc.harmonic.r,xi_max, ...
        'opposing','tuned');
    G = sys.Gamma(1);

    expirTau = exp(1i*linspace(0,2*pi,sol.N_Sample));
    
    Qloc = real(Qana(1)*expirTau);
    Uloc = real(1i*exc.harmonic.r*Qana(1)*expirTau);

    Qnonloc = real(Qana(4)*expirTau);
    Unonloc = real(1i*exc.harmonic.r*Qana(4)*expirTau);

    figure(5)
    stem(0:(sys.N_s-1),abs(Qana)/sys.qref,'-+','Color',color.analytics,...
        'LineWidth',1,'MarkerSize',15,'DisplayName','Analytical')
    hold on;
    plot(0:(sys.N_s-1),qhat/sys.qref,'-o','Color',color.ies,...
        'LineWidth',1.5,'DisplayName','Num. Sim.','MarkerSize',8)
    box on;
    plot([0,sys.N_s/2],qhat([1,sys.N_s/2+1])/sys.qref,'o','Color',color.ies,...
        'LineWidth',1.5,'MarkerFaceColor',color.ies,'MarkerSize',8,...
        'HandleVisibility','off')
    xlabel('Sector - $j$')
    ylabel('$\hat{q}_j/\hat{q}_\mathrm{ref}$')
    xlim([-0.499 (sys.N_s-1)+0.499])
    ylim([0.95*min(qhat), 1.05*max(qhat)]/sys.qref)
    legend;
    savefig([savepath 'amplitude_distribution_example.fig'])


    
    figure(6)
    p1 = plot(Q(1,:)/sys.qref,U(1,:)/sys.qref/exc.harmonic.r,...
        'LineWidth',0.5,'Color',color.ies);
    hold on;
    p1.Color(4) = 0.05;
    plot(Qloc/sys.qref,Uloc/sys.qref/exc.harmonic.r,...
        '--','LineWidth',2.5,'Color',color.analytics);
    box on;
    title('Synchronized Sector')
    xlabel('$q/\hat{q}_\mathrm{ref}$')
    ylabel('$u/\hat{u}_\mathrm{ref}$')
    savefig([savepath 'phasespace_synch_example.fig'])

    figure(7)
    p2 = plot(Q(4,:)/sys.qref,U(4,:)/sys.qref/exc.harmonic.r,...
        'LineWidth',0.5,'Color',color.ies);
    hold on;
    p2.Color(4) = 0.05;
    plot(Qnonloc/sys.qref,Unonloc/sys.qref/exc.harmonic.r,...
        '--','LineWidth',2.5,'Color',color.analytics);
    box on;
    title('Non-Synchronized Sector')
    xlabel('$q/\hat{q}_\mathrm{ref}$')
    ylabel('$u/\hat{u}_\mathrm{ref}$')
    savefig([savepath 'phasespace_nonsynch_tuned_example.fig'])

    % Plot absorber motion
    figure(8)
    plot(exc.harmonic.r*TAU/2/pi,QA(1,:),'LineWidth',1,...
        'Color',color.iesabsorber)
    hold on;
    plot(exc.harmonic.r*TAU/2/pi,Q(1,:)+sys.Gamma(1),'LineWidth',1,...
        'Color',color.ies)
    plot(exc.harmonic.r*TAU/2/pi,Q(1,:)-sys.Gamma(1),'LineWidth',1,...
        'Color',color.ies)
    box on;
    xlabel('$r\tau / (2 \pi)$')
    ylabel('$q_\mathrm{a}$')
    axis tight;
    title('Synchronized Sector')
    savefig([savepath 'absmove_synch_tuned_example.fig'])

    figure(9)
    plot(exc.harmonic.r*TAU/2/pi,QA(4,:),'LineWidth',1,...
        'Color',color.iesabsorber)
    hold on;
    plot(exc.harmonic.r*TAU/2/pi,Q(4,:)+sys.Gamma(4),'LineWidth',1,...
        'Color',color.ies)
    plot(exc.harmonic.r*TAU/2/pi,Q(4,:)-sys.Gamma(4),'LineWidth',1,...
        'Color',color.ies)
    box on;
    xlabel('$r\tau / (2 \pi)$')
    ylabel('$q_\mathrm{a}$')
    axis tight;
    title('Non-Synchronized Sector')
    savefig([savepath 'absmove_nonsynch_tuned_example.fig'])
end