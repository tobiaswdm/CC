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

%% Compute tuned manifold

% Auxilliary Variable
rho = (2/pi) * (1-sys.eN) / (1+sys.eN);
% Minimum amplitude - turning point of SIM
% Including safety of 0.1% higher amplitude
xi_min = 1.001 * rho/sqrt(1+rho^2);

% Clearance normalized amplitudes
xi = logspace(log10(xi_min),...
    log10(simsetup.SynchronizationSingleSectorStability.xi_max),...
    simsetup.SynchronizationSingleSectorStability.Nxi);
r = linspace(simsetup.SynchronizationSingleSectorStability.r_range(1),...
    simsetup.SynchronizationSingleSectorStability.r_range(2), ...
    simsetup.SynchronizationSingleSectorStability.Nr);

% Linear FRF
q_fixed = abs(ComputeLinearResponse(r,sys,exc,'tuned','fixed_absorbers'));
q_fixed = q_fixed(1,:);
q_removed = abs(ComputeLinearResponse(r,sys,exc,'tuned','removed_absorbers'));
q_removed = q_removed(1,:);

% FRS
[Gamma_Scale,Xi,R] = SingleSectorFRS(xi,r,sys,exc,'tuned');

% Plot FRS
figure(2);
surf(R,Xi,Gamma_Scale,'EdgeColor','none')
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
c = contourc(r,xi,Gamma_Scale',[sys.Gamma_Scale sys.Gamma_Scale]);

% Determine max amplitude FRF
[qhat_max,qhat_max_violated,qhat_syn,r_plot] = ...
    LocalizedFrequencyAmplitudeCurve(c,sys,exc,'single','tuned');

% Coarsen contour for stability analysis
c = CoarsenContour(c,...
    simsetup.SynchronizationSingleSectorStability.stepsize);

% Study asymptotic and practical stability of tuned system
[qhat_practically_stable_t,qhat_stable_t,qhat_unstable_t, ... 
qhatsynch_practically_stable_t,qhatsynch_stable_t,qhatsynch_unstable_t, ...
IPR_t,LF_t,r_num_t] = StabilityAnalysisLSR(c,sys,sol,exc, ...
                                            'single','tuned',true,true);

% Get maximum practically stable amplitude in tuned case
if ~isempty(qhat_practically_stable_t)
    qhat_practically_stable_max_t = max(qhat_practically_stable_t);
else
    qhat_practically_stable_max_t = NaN;
end


figure(4);
hold on;
plot(r,q_fixed/sys.qref,...
    'LineWidth',.5,'Color',color.reference,'DisplayName', ...
    'Fixed VI-NESs T')
plot(r,q_removed/sys.qref,'-.',...
    'LineWidth',.5,'Color',color.reference,'DisplayName', ...
    'Removed VI-NESs T')
if simsetup.SynchronizationSingleSectorStability.N_MCS == 0
    plot(r_plot,qhat_syn/sys.qref,'--',...
            'LineWidth',1.5,'Color',color.background,'DisplayName', ...
            'Syn. sector')
    plot(r_plot,qhat_max/sys.qref,...
            'LineWidth',1.5,'Color',color.ies,'DisplayName', ...
            'Max. amp.')
    plot(r_plot,qhat_max_violated/sys.qref,':',...
                'LineWidth',1.5,'Color',color.show,'DisplayName', ...
                'Viol. kin. constr.')
else
    plot(r_plot,qhat_max/sys.qref,...
            'LineWidth',1.5,'Color',color.ies,'DisplayName', ...
            'LSR Analytical T')
    plot(r_plot,qhat_max_violated/sys.qref,':',...
            'LineWidth',1.5,'Color',color.show,'DisplayName', ...
            'Viol. kin. constr. T')
end
set(gca,'YScale','log')
axis tight;
box on;
xlabel('$r$')
ylabel('$\hat{q}/\hat{q}_\mathrm{ref}$')


%% Study realizations of mistuned systems

% Count if practically stable solutions were found
k = 0;

% Save maximum practically stable amplitude
qhat_practically_stable_max_mt = nan(1, ...
    simsetup.SynchronizationSingleSectorStability.N_MCS);
% Save IPR of maximum practically stable amplitude
IPR_practically_stable_max_mt = nan(1, ...
    simsetup.SynchronizationSingleSectorStability.N_MCS);
% Save LF of maximum practically stable amplitude
LF_practically_stable_max_mt = nan(1, ...
    simsetup.SynchronizationSingleSectorStability.N_MCS);
if sys.sigma_omega ~= 0
% Save corresponding frequency of maximum practically stable amplitude
    r_practically_stable_max_mt = nan(1, ...
        simsetup.SynchronizationSingleSectorStability.N_MCS);
% Save local eigenfrequency mistuning of sector with largest practically
% stable amplitude during LSR
    delta_omega_practically_stable_max = nan(1, ...
        simsetup.SynchronizationSingleSectorStability.N_MCS);
% Check if sector with largest practically stable amplitude during LSR also
% has the lowest local eigenfrequency for the mistuning pattern
    isminfreq_practically_stable_max = nan(1, ...
        simsetup.SynchronizationSingleSectorStability.N_MCS);
end

% Cell Arrays for amplitudes, frequencies and localization measures
qhat_practically_stable = cell(1, ...
    simsetup.SynchronizationSingleSectorStability.N_MCS);
r_num = cell(1, simsetup.SynchronizationSingleSectorStability.N_MCS);
LF = cell(1, simsetup.SynchronizationSingleSectorStability.N_MCS);
IPR = cell(1, simsetup.SynchronizationSingleSectorStability.N_MCS);

% Store if specific LSR has practically stable solution
ispracticallystable = false(sys.N_s, ...
    simsetup.SynchronizationSingleSectorStability.N_MCS);

for i = 1:simsetup.SynchronizationSingleSectorStability.N_MCS

    disp(['Testing mistuned system ' num2str(i) ' of '...
        num2str(simsetup.SynchronizationSingleSectorStability.N_MCS)])
    
    % Intialize cell array for practically stable amplitudes and freqs.
    qhat_practically_stable_array = cell(1,sys.N_s);
    r_num_array = cell(1,sys.N_s);
    % Intialize cell array for localization measures
    LF_array = cell(1,sys.N_s);
    IPR_array = cell(1,sys.N_s);
    % Maximum practically stable amplitude from circular shifts
    qhat_practically_stable_max = zeros(1,sys.N_s);
    r_practically_stable_max = zeros(1,sys.N_s);
    % Localization measures of maximum amplitude
    IPR_practically_stable_max = zeros(1,sys.N_s);
    LF_practically_stable_max = zeros(1,sys.N_s);
    
    % Loop through circular shifts of mistuning patterns
    % This is equivalent as shifting the localized sector j=0,...,(Ns-1)
    for j = 1:sys.N_s
        
        if j == 1
            % Build new mistuned system
            [sys_mt,exc_mt] = BuildSystem(sys,exc,'mistuned');
        else
            % Circular shift of mistuning patterns by one sector
            sys_mt.delta_omega = circshift(sys_mt.delta_omega,1,1);
            sys_mt.delta_g = circshift(sys_mt.delta_g,1,1);

            % Rebuild mistuned system
            [sys_mt,exc_mt] = BuildSystem(sys_mt,exc,'mistuned_defined');
        end

        % Determine mistuned FRS
        [Gamma_Scale_mt,~,~] = ...
            SingleSectorFRS(xi,r,sys_mt,exc_mt,'mistuned');

        % Get Level curves at clearance
        c = contourc(r,xi,Gamma_Scale_mt',...
            [sys_mt.Gamma_Scale sys_mt.Gamma_Scale]*(1+sys_mt.delta_g(1)));
    
        % Coarsen contour for stability analysis
        c = CoarsenContour(c,...
            simsetup.SynchronizationSingleSectorStability.stepsize);
    
        % Study only practical stability
        [qhat_practically_stable_array{j},~,~,~,~,~,IPR_array{j},LF_array{j},...
            r_num_array{j}] =StabilityAnalysisLSR(c,sys_mt,sol,exc_mt, ...
            'single','mistuned',true,true);
        
        % Extract maximum
        if isempty(qhat_practically_stable_array{j})
            % No solution found
            qhat_practically_stable_max(j) = NaN;
            r_practically_stable_max(j) = NaN;
            IPR_practically_stable_max(j) = NaN;
            LF_practically_stable_max(j) = NaN;
        else
            % Extract maximum
            [qhat_practically_stable_max(j),imax] = ...
                max(qhat_practically_stable_array{j});
            r_practically_stable_max(j) = r_num_array{j}(imax);
            IPR_practically_stable_max(j) = IPR_array{j}(imax);
            LF_practically_stable_max(j) = LF_array{j}(imax);

            % Save that practically stable solution was found
            ispracticallystable(j,i) = true;
        end
        

    end
    
    % Practically stable solution found?
    if any(~isnan(qhat_practically_stable_max))

        % Apply circular shift to obtain initial configuration
        sys_mt.delta_omega = circshift(sys_mt.delta_omega,1,1);

        % Increase count
        k = k+1;

        % Extract shift config with maximum pract. stable amplitude and
        % corresponding measures
        [qhat_practically_stable_max_mt(i),j_max] = ...
        max(qhat_practically_stable_max);
        if sys_mt.sigma_omega ~= 0
            r_practically_stable_max_mt(i) = ...
                    r_practically_stable_max(j_max);
            delta_omega_practically_stable_max(i) = ...
                    sys_mt.delta_omega(j_max);
            isminfreq_practically_stable_max(i) = ...
                    sys_mt.delta_omega(j_max) == ...
                    min(sys_mt.delta_omega);
            IPR_practically_stable_max_mt(i) = ...
                    IPR_practically_stable_max(j_max);
            LF_practically_stable_max_mt(i) = ...
                    LF_practically_stable_max(j_max);
        end
    else
        j_max  = 1;
        qhat_practically_stable_max_mt(i) = NaN;
        if sys_mt.sigma_omega ~= 0
            r_practically_stable_max_mt(i) = NaN;
            delta_omega_practically_stable_max(i) = NaN;
            isminfreq_practically_stable_max(i) = NaN;
            IPR_practically_stable_max_mt(i) = NaN;
            LF_practically_stable_max_mt(i) = NaN;
        end
    end
    
    % Pick largest practically stable amplitudes
    qhat_practically_stable{i} = qhat_practically_stable_array{j_max};
    r_num{i} = r_num_array{j_max};
    LF{i} = LF_array{j_max};
    IPR{i} = IPR_array{j_max};
end

save([savepath 'r_num.mat'],'r_num')
save([savepath 'qhat_practically_stable.mat'], ...
    'qhat_practically_stable')
save([savepath 'qhat_practically_stable_max_mt.mat'], ...
    'qhat_practically_stable_max_mt')
save([savepath 'IPR_practically_stable_max_mt.mat'], ...
    'IPR_practically_stable_max_mt')
save([savepath 'LF_practically_stable_max_mt.mat'], ...
    'LF_practically_stable_max_mt')
save([savepath 'IPR.mat'],'IPR')
save([savepath 'LF.mat'],'LF')

% Convert to matrices
qhat_practically_stable = cell2mat(qhat_practically_stable);
r_num = cell2mat(r_num);
switch simsetup.SynchronizationSingleSectorStability.LocalizationMeasure
    case 'LF'
        locmeasure = cell2mat(LF);
    case 'IPR'
        locmeasure = cell2mat(IPR);
end


if simsetup.SynchronizationSingleSectorStability.N_MCS>0
    figure(4);
    hold on;
    scatter(r_num,qhat_practically_stable/sys.qref,20,locmeasure,'filled',...
    'MarkerFaceAlpha',1,'MarkerEdgeAlpha',0,...
    'Displayname','Pract. Stable MT','Marker','square')
    colormap(flipud(gray).^0.8);
    cb = colorbar();
    switch simsetup.SynchronizationSingleSectorStability.LocalizationMeasure
        case 'LF'
            ylabel(cb,'LF','Rotation',90,'Interpreter','latex')
        case 'IPR'
            ylabel(cb,'IPR','Rotation',90,'Interpreter','latex')
    end
end



% In how many mistuned cases were practically stable solutions found?
prac_stab_ratio = k/simsetup.SynchronizationSingleSectorStability.N_MCS;
save([savepath 'prac_stab_ratio.mat'],'prac_stab_ratio')
disp(['Practically stable solutions found in ' ...
    num2str(100*prac_stab_ratio) ' percent of the mistuned cases.'])

% In how many mistuned cases were practically stable solutions found under
% consdieration of different synchronized sectors for same mistuning pattern?
prac_stab_ratio_total = sum(ispracticallystable,'all')/...
                  numel(ispracticallystable);
save([savepath 'prac_stab_prac_stab_ratio_total.mat'],'prac_stab_ratio_total')
disp(['Practically stable solutions found in ' ...
    num2str(100*prac_stab_ratio_total) ' percent of the mistuned cases and' ...
    ' synchronized sectors.'])

figure(4)
hold on;
if simsetup.SynchronizationSingleSectorStability.N_MCS == 0
    scatter(r_num_t,qhatsynch_unstable_t/sys.qref,20,'MarkerFaceColor',color.show,...
    'MarkerEdgeColor','k','Displayname','Unstable')
    scatter(r_num_t,qhatsynch_stable_t/sys.qref,20,'MarkerFaceColor',myColors('cyan'),...
        'MarkerEdgeColor','k','Displayname','Stable')
    scatter(r_num_t,qhatsynch_practically_stable_t/sys.qref,40,'pentagram',...
        'MarkerFaceColor',myColors('green'),'MarkerEdgeColor','k',...
        'Displayname','Pract. Stable')
else
    scatter(r_num_t,qhat_unstable_t/sys.qref,20,'MarkerFaceColor',color.show,...
    'MarkerEdgeColor','k','Displayname','Unstable T')
    scatter(r_num_t,qhat_stable_t/sys.qref,20,'MarkerFaceColor',myColors('cyan'),...
        'MarkerEdgeColor','k','Displayname','Stable T')
    scatter(r_num_t,qhat_practically_stable_t/sys.qref,40,'pentagram',...
        'MarkerFaceColor',myColors('green'),'MarkerEdgeColor','k',...
        'Displayname','Pract. Stable T')
end
legend;
savefig([savepath 'frequency_amplitude_stability.fig'])

%% Show example of full mistuned FRF and stability

if (sys.sigma_g ~= 0 || sys.sigma_omega ~= 0) || ...
    (isfield(sys,'delta_g') || isfield(sys,'delta__omega'))    

    % Build mistuned system
    if isfield(sys,'delta_g') && isfield(sys,'delta_omega')
        [sys_mt,exc_mt] = BuildSystem(sys,exc,'mistuned_defined');
    else
        [sys_mt,exc_mt] = BuildSystem(sys,exc,'mistuned');
    end
    % Determine mistuned FRS
    [Gamma_Scale_mt,~,~] = ...
        SingleSectorFRS(xi,r,sys_mt,exc_mt,'mistuned');

    % Save mistuning configuration
    delta_omega = sys_mt.delta_omega;
    delta_g = sys_mt.delta_g;
    save([savepath 'delta_omega_example.mat'], 'delta_omega')
    save([savepath 'delta_g_example.mat'], 'delta_g')
    
    % Get Level curves at clearance
    c = contourc(r,xi,Gamma_Scale_mt',...
        [sys_mt.Gamma_Scale sys_mt.Gamma_Scale]*(1+sys_mt.delta_g(1)));
    
    % Determine max amplitude FRF
    [qhat_max,qhat_max_violated,qhat_localized,r_plot] = ...
        LocalizedFrequencyAmplitudeCurve(c,sys_mt,exc_mt,'single','mistuned');
    
    % Coarsen contour for stability analysis
    c = CoarsenContour(c,...
        simsetup.SynchronizationSingleSectorStability.stepsize);
    
    % Study asymptotic and practical stability of mistuned system
    [qhat_practically_stable,qhat_stable,qhat_unstable,~,~,~,~,~,r_num] =...
        StabilityAnalysisLSR(c,sys_mt,sol,exc_mt,'mistuned',true,true);
    
    % Plot FRF
    figure(5);
    hold on;
    plot(r,q_fixed/sys.qref,...
        'LineWidth',.5,'Color',color.reference,'DisplayName', ...
        'Fixed abs.')
    plot(r,q_removed/sys.qref,'-.',...
        'LineWidth',.5,'Color',color.reference,'DisplayName', ...
        'Removed abs.')
    plot(r_plot,qhat_localized/sys.qref,'--',...
        'LineWidth',1.5,'Color',color.background,'DisplayName', ...
        'Synch. sector')
    plot(r_plot,qhat_max/sys.qref,...
                'LineWidth',1.5,'Color',color.ies,'DisplayName', ...
                'Tuned')
    plot(r_plot,qhat_max_violated/sys.qref,':',...
                'LineWidth',1.5,'Color',color.show,'DisplayName', ...
                'Tuned - Viol. kin. constr.')
    scatter(r_num,qhat_unstable/sys.qref,20,'MarkerFaceColor',color.show,...
        'MarkerEdgeColor','k','Displayname','Unstable')
    scatter(r_num,qhat_stable/sys.qref,20,'MarkerFaceColor',myColors('cyan'),...
        'MarkerEdgeColor','k','Displayname','Stable')
    scatter(r_num,qhat_practically_stable/sys.qref,40,'pentagram',...
        'MarkerFaceColor',myColors('green'),'MarkerEdgeColor','k',...
        'Displayname','Pract. Stable')
    set(gca,'YScale','log')
    title('Mistuned System')
    axis tight;
    legend;
    box on;
    xlabel('$r$')
    ylabel('$\hat{q}/\hat{q}_\mathrm{ref}$')
    savefig([savepath 'frequency_amplitude_mistuned_example.fig'])
    
    % Plot FRS of mistuned system
    figure(6);
    surf(R,Xi,Gamma_Scale_mt,'EdgeColor','none','DisplayName','FRS')
    hold on;
    contour3(R,Xi,Gamma_Scale_mt,[1 1]*sys.Gamma_Scale, '-k', 'LineWidth',3,...
        'DisplayName','Nominal Clearance')
    if sys.sigma_g ~= 0 % Only if clearance is mistuned
       contour3(R,Xi,Gamma_Scale_mt,[1 1]*(1+sys.sigma_g)*sys.Gamma_Scale, ...
           '--k', 'LineWidth',1.5,'DisplayName','Mistuned Clearance Std.')
       contour3(R,Xi,Gamma_Scale_mt,[1 1]*(1-sys.sigma_g)*sys.Gamma_Scale, ...
           '--k', 'LineWidth',1.5,'HandleVisibility','off')
    end
    title('Mistuned System')
    box on;
    colormap turbo
    xlabel('$r$')
    ylabel('$\xi$')
    zlabel('$\Gamma/\hat{q}_\mathrm{ref}$')
    set(gca,'YScale','log')
    axis tight;
end

%% Highest practically stable amplitude in tuned system

if any(~isnan(qhat_practically_stable_t))
    
    % Extract highest practically stable amplitude
    %[xi_max_t,i_max_t] = max(qhat_practically_stable_t/sys.Gamma(1));
    % Extract highest practically stable amplitude
    [xi_max_t,i_max_t] = min(qhat_practically_stable_t/sys.Gamma(1));
    exc.harmonic.r = r_num_t(i_max_t);

    % Get initial conditions for time simulation
    % Get initial conditions
    sol.N_Tau = 300;
    [sol.q0,sol.u0,sol.qa0,sol.ua0] = ...
        LocalizedInitialConditions(sys,exc,xi_max_t,exc.harmonic.r,...
        'single','tuned','practical_stability');
    sol = ConfigureIntegrator(sol,sys,exc,...
                    'no change',false,'tuned');
    sol.NP_trans = sol.N_P*1000;

    % Simulation
    [ETA,QA,Chi,~,TAU] = MoreauIntegration(sys,exc,sol,'tuned');
    Q = sys.Phi*ETA;
    U = sys.Phi*Chi;

    % Extract amplitudes
    [qhat,~] = MeanAmplitude(Q(:,2:end),sol.N_Sample);
    
    figure(7)
    plot(0:(sys.N_s-1),qhat/sys.qref,'-o','Color',color.ies,...
        'LineWidth',1.5,'DisplayName','Num. Sim.','MarkerSize',8)
    hold on;
    box on;
    plot(0,qhat(1)/sys.qref,'o','Color',color.ies,...
        'LineWidth',1.5,'MarkerFaceColor',color.ies,'MarkerSize',8,...
        'HandleVisibility','off')
    xlabel('Sector - $j$')
    ylabel('$\hat{q}_j/\hat{q}_\mathrm{ref}$')
    xlim([-0.499 (sys.N_s-1)+0.499])
    ylim([0.95*min(qhat), 1.05*max(qhat)]/sys.qref)


    % Average modal energies tuned system
    [E_mod,E_mod_avg] = ModalEnergies(sys,sol,ETA,Chi);

    figure(8)
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
    savefig([savepath 'modal energies_tuned_example.fig'])
    
    % Phase space tuned system
    % Get analytical HB approximation
    [Qana,~,~] = ...
        RecoverCondensedDOFs(sys,exc,exc.harmonic.r,xi_max_t, ...
        'single','tuned');
    G = sys.Gamma(1);

    expirTau = exp(1i*linspace(0,2*pi,sol.N_Sample));
    
    Qloc = real(Qana(1)*expirTau);
    Uloc = real(1i*exc.harmonic.r*Qana(1)*expirTau);

    Qnonloc = real(Qana(4)*expirTau);
    Unonloc = real(1i*exc.harmonic.r*Qana(4)*expirTau);

    figure(7)
    stem(0:(sys.N_s-1),abs(Qana)/sys.qref,'-+','Color',color.analytics,...
        'LineWidth',1,'MarkerSize',15,'DisplayName','Analytical')
    hold on;
    plot(0:(sys.N_s-1),qhat/sys.qref,'-o','Color',color.ies,...
        'LineWidth',1.5,'DisplayName','Num. Sim.','MarkerSize',8)
    box on;
    plot(0,qhat(1)/sys.qref,'o','Color',color.ies,...
        'LineWidth',1.5,'MarkerFaceColor',color.ies,'MarkerSize',8,...
        'HandleVisibility','off')
    xlabel('Sector - $j$')
    ylabel('$\hat{q}_j/\hat{q}_\mathrm{ref}$')
    xlim([-0.499 (sys.N_s-1)+0.499])
    ylim([0.95*min(qhat), 1.05*max(qhat)]/sys.qref)
    legend;
    savefig([savepath 'amplitude_distribution_tuned_example.fig'])

    
    figure(9)
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
    savefig([savepath 'phasespace_synch_tuned_example.fig'])

    figure(10)
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
    figure(11)
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

    figure(12)
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

% Plot PDF of maximum amplitude
if any(~isnan(qhat_practically_stable_max_mt)) && ...
        simsetup.SynchronizationSingleSectorStability.N_MCS~=0
    
    [E,x] = ecdf(qhat_practically_stable_max_mt/sys.qref);
    
    figure(13)
    hold on;
    colororder({color.ies;color.analytics})
    yyaxis right;
    stairs(x,E,'--k','LineWidth',1.5)
    ylabel('CDF')
    yyaxis left;
    histogram(qhat_practically_stable_max_mt/sys.qref,...
    'Normalization','pdf','FaceColor',color.ies,'LineStyle','none');
    box on;
    xlabel('$\mathrm{max}_r \left\{ \hat{q}^\ast / \hat{q}_\mathrm{ref} \right\}$')
    ylabel('PDF')
    axis tight
    savefig([savepath 'qhat_max_pdf.fig'])
    
     % Determine performance of tuned system with absorber
    
    % Frequency range for stepping
    r_range = [sys.r_k(exc.k+1), sys.r_k_noabs(exc.k+1)];
    r_range = simsetup.SynchronizationSingleSectorStability.r_scale.*r_range;
    r_steps = linspace(r_range(1),r_range(2),...
                       simsetup.SynchronizationSingleSectorStability.N_rSteps);
    
    % Resonance amplitude tuned system
    [qhat_res_tuned, ~, ~] = ....
    FindResonance(sys,sol,exc,r_steps,'tuned');
    
    [E,x] = ecdf(qhat_practically_stable_max_mt/qhat_res_tuned);
    
    figure(14)
    hold on;
    colororder({color.ies;color.analytics})
    yyaxis right;
    stairs(x,E,'--k','LineWidth',1.5)
    ylabel('CDF')
    yyaxis left;
    histogram(qhat_practically_stable_max_mt/qhat_res_tuned,...
    'Normalization','pdf','FaceColor',color.ies,'LineStyle','none');
    box on;
    xlabel('$A$')
    ylabel('PDF')
    axis tight
    savefig([savepath 'A_max_pdf.fig'])
end

% Plot correlation between frequency of highest practically stable
% amplitude and local eigenfrequency mistuning
if sys.sigma_omega ~= 0
    
    save([savepath 'r_practically_stable_max_mt.mat'], ...
        'r_practically_stable_max_mt')
    save([savepath 'delta_omega_practically_stable_max.mat'], ...
        'delta_omega_practically_stable_max')
    save([savepath 'isminfreq_practically_stable_max.mat'], ...
        'isminfreq_practically_stable_max')

    figure(15);
    hold on; box on;
    switch simsetup.SynchronizationSingleSectorStability.LocalizationMeasure
        case 'LF'
            scatter(r_practically_stable_max_mt, ...
            delta_omega_practically_stable_max, ...
            20,LF_practically_stable_max_mt,'filled',...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',0,...
            'Displayname','Pract. Stable MT','Marker','square')
            colormap(flipud(gray).^0.8);
            cb = colorbar();
            ylabel(cb,'LF','Rotation',90,'Interpreter','latex')
        case 'IPR'
            scatter(r_practically_stable_max_mt, ...
            delta_omega_practically_stable_max, ...
            20,IPR_practically_stable_max_mt,'filled',...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',0,...
            'Displayname','Pract. Stable MT','Marker','square')
            colormap(flipud(gray).^0.8);
            cb = colorbar();
            ylabel(cb,'IPR','Rotation',90,'Interpreter','latex')
    end
    xlabel('$r_\mathrm{max}$')
    ylabel('$\delta_{r,j}$')
    title(['PCC = ' num2str(round(1000*corr( ...
        r_practically_stable_max_mt(~isnan( ...
        delta_omega_practically_stable_max))',...
        delta_omega_practically_stable_max(~isnan( ...
        delta_omega_practically_stable_max))'))/1000)], ...
        'Interpreter','none')
    savefig([savepath 'r_practically_stable_max_mt.fig'])
    
    figure(16);
    hold on; box on;
    histogram(isminfreq_practically_stable_max,'FaceColor',color.ies, ...
        'EdgeColor','k','Normalization','probability')
    title('Highest amplitude in softest sector?')
    ylabel('Probability')
    xticks([0 1])
    xticklabels({'No','Yes'})
    axis tight;
    ylim([0 1])
    savefig([savepath 'isminfreq_practically_stable_max.fig'])

end
