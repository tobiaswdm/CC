%% Build the system

% Tuned
[sys,exc] = BuildSystem(sys,exc,'tuned');

% Draw Dispersion Diagram
DrawDispersion(sys,color,savepath);

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
[ETA_mt,QA_mt,Chi_mt,UA_mt,~] = MoreauIntegration(sys,exc,sol_mt,'mistuned');
Q_mt = sys.Phi * ETA_mt;
U_mt = sys.Phi * Chi_mt;

% Simulation tuned system
[ETA,QA,Chi,UA,TAU] = MoreauIntegration(sys,exc,sol,'tuned');
Q = sys.Phi * ETA;
U = sys.Phi * Chi;
TAU = TAU - TAU(1);

%% Post-processing

% Amplitudes - mistuned system
[qhat_mt,qhat_mt_std] = MeanAmplitude(Q_mt(:,2:end),sol_mt.N_Sample);

% Amplitudes - tuned system
[qhat,qhat_std] = MeanAmplitude(Q(:,2:end),sol.N_Sample);

% Modal energies
[E_mod] = ModalEnegies(sys,sol,ETA,Chi);
[E_mod_mt] = ModalEnegies(sys,sol_mt,ETA_mt,Chi_mt);

% Spatial energies
E = SpatialEnergies(sys,sol,Q,U,UA,'tuned');
E_mt = SpatialEnergies(sys,sol,Q_mt,U_mt,UA_mt,'mistuned');

% Hilbert transform
QH = hilbert(Q')';
Theta = unwrap(angle(QH),[],2);
Theta = diff([Theta;Theta(1,:)],[],1);

%% Figures

figure(2);
stem(0:(sys.N_s-1),qhat/sys.qref,"filled",'Color',color.background, ...
    'LineWidth',2,'MarkerSize',10,'DisplayName','Tuned System')
hold on;
errorbar(0:(sys.N_s-1),qhat/sys.qref,...
    qhat_std/sys.qref,"LineStyle","none",...
    'LineWidth',1.5,'CapSize',6,...
    'Color',color.background,...
    'HandleVisibility','off')
stem(0:(sys.N_s-1),qhat_mt/sys.qref,"filled",'Color',color.ies, ...
    'LineWidth',2,'MarkerSize',10,'DisplayName','Mistuned System')
errorbar(0:(sys.N_s-1),qhat_mt/sys.qref,...
    qhat_mt_std/sys.qref,"LineStyle","none",...
    'LineWidth',1.5,'CapSize',10,...
    'Color',color.ies,...
    'HandleVisibility','off')
xlabel('Sector - $j$')
ylabel('$\hat{q}_j/\hat{q}_\mathrm{ref}$')
yline(mean(qhat)/sys.qref,'--k','DisplayName','$\hat{q}_\mathrm{mean}$')
legend;
xlim([-0.5 sys.N_s-0.5])
savefig([savepath 'amplitude_distribution.fig'])

% Tuned Displacement
index=[1 6 2 7 3 8 4 9 5 10];
figure(3);
hold on;
title('Displacement')
for i = 1:min(sys.N_s,10)
    name = ['$q_' num2str(index(i)-1) '$'];
    subplot(5,2,i)
    hold on;
    box on;
    plot(exc.harmonic.r*TAU/2/pi,Q(index(i),:),'LineWidth',1,'Color',color.ies)
    yline(qhat(index(i)),'-','LineWidth',.5,'Color',color.show)
    yline(qhat(index(i))+qhat_std(index(i)),'--','LineWidth',.5,'Color',color.show)
    xlim(exc.harmonic.r*[TAU(1) TAU(end)]/2/pi)
    %xlim([TAU(end)-30*(2*pi/exc.harmonic.r) TAU(end)])
    xlabel('$r\tau / (2 \pi)$')
    ylim([-1.3 1.3]*max(qhat+2*qhat_std))
    ylabel(name)
end

% Mistuned Displacement
figure(4);
hold on;
title('Displacement - Mistuned')
for i = 1:min(sys.N_s,10)
    name = ['$q_' num2str(index(i)-1) '^\ast $'];
    subplot(5,2,i)
    hold on;
    box on;
    plot(exc.harmonic.r*TAU/2/pi,Q_mt(index(i),:),'LineWidth',1,'Color',color.ies)
    yline(qhat_mt(index(i)),'-','LineWidth',.5,'Color',color.show)
    yline(qhat_mt(index(i))+qhat_mt_std(index(i)),'--','LineWidth',.5,'Color',color.show)
    xlim(exc.harmonic.r*[TAU(1) TAU(end)]/2/pi)
    %xlim([TAU(end)-30*(2*pi/exc.harmonic.r) TAU(end)])
    xlabel('$r\tau / (2 \pi)$')
    ylim([-1.3 1.3]*max(qhat_mt+2*qhat_mt_std))
    ylabel(name)
end

% Absorber movement tuned
figure(5);
hold on;
for i = 1:min(sys.N_s,10)
    name = ['$q_{\mathrm{a},' num2str(index(i)-1) '}$'];
    subplot(5,2,i)
    hold on;
    box on;
    plot(exc.harmonic.r*TAU/2/pi,...
        Q(index(i),:)+sys.Gamma(2*index(i)),...
        'LineWidth',1,'Color',color.ies)
    plot(exc.harmonic.r*TAU/2/pi,...
        Q(index(i),:)-sys.Gamma(2*index(i)),...
        'LineWidth',1,'Color',color.ies)
    plot(exc.harmonic.r*TAU/2/pi,...
        QA(index(i),:),'LineWidth',1,'Color',color.iesabsorber)
    xlim(exc.harmonic.r*[TAU(1) TAU(end)]/2/pi)
    %xlim([TAU(end)-30*(2*pi/exc.harmonic.r) TAU(end)])
    xlabel('$r\tau / (2 \pi)$')
    ylim([-1.3 1.3]*(max(qhat_mt+2*qhat_mt_std+sys.Gamma(2*index(i)))))
    ylabel(name)
end

% Absorber movement mistuned
figure(6);
hold on;
for i = 1:min(sys.N_s,10)
    name = ['$q_{\mathrm{a},' num2str(index(i)-1) '}^\ast $'];
    subplot(5,2,i)
    hold on;
    box on;
    plot(exc.harmonic.r*TAU/2/pi,...
        Q_mt(index(i),:)+sys.Gamma_mt(2*index(i)),...
        'LineWidth',1,'Color',color.ies)
    plot(exc.harmonic.r*TAU/2/pi,...
        Q_mt(index(i),:)-sys.Gamma_mt(2*index(i)),...
        'LineWidth',1,'Color',color.ies)
    plot(exc.harmonic.r*TAU/2/pi,...
        QA_mt(index(i),:),'LineWidth',1,'Color',color.iesabsorber)
    xlim(exc.harmonic.r*[TAU(1) TAU(end)]/2/pi)
    %xlim([TAU(end)-30*(2*pi/exc.harmonic.r) TAU(end)])
    xlabel('$r\tau / (2 \pi)$')
    ylim([-1.3 1.3]*(max(qhat_mt+2*qhat_mt_std+max(sys.Gamma_mt))))
    ylabel(name)
end

% Averaged Sector Energies Tuned
figure(7);
imagesc(exc.harmonic.r*TAU/2/pi,0:(sys.N_s-1),E./sum(E,1))
hold on;
colormap(1-pink.^2)
cb = colorbar(); 
ylabel(cb,'$E_j^\ast / E (\tau) $','Rotation',90,'Interpreter','latex')
title('Averaged Sector Energies Tuned')
xlabel('$r\tau / (2 \pi)$')
ylabel('Sector $j$')
ylabel('$E_j (\tau) / E (\tau) $')

% Averaged Sector Energies Mistuned
figure(8)
imagesc(exc.harmonic.r*TAU/2/pi,0:(sys.N_s-1),E_mt)
hold on;
colormap(1-pink.^2)
cb = colorbar(); 
ylabel(cb,'$E_j^\ast / E (\tau) $','Rotation',90,'Interpreter','latex')
title('Averaged Sector Energies Mistuned')
xlabel('$r\tau / (2 \pi)$')
ylabel('Sector $j$')

% State space tuned system
figure(9)
plot3(Q(1,:),Q(2,:),Q(3,:),'DisplayName','$j=0$',...
    'LineWidth',0.5,'Color',color.ies);
hold on;
box on;
axis tight;
xlabel('$q_j$')
ylabel('$q_{j+1}$')
zlabel('$q_{j+2}$')
grid on;
legend;
title('State Space Tuned system')

% State space mistuned system
figure(10)
plot3(Q_mt(1,:),Q_mt(2,:),Q_mt(3,:),'DisplayName','$j=0$',...
    'LineWidth',0.5,'Color',color.ies);
hold on;
box on;
axis tight;
xlabel('$q^\ast_j$')
ylabel('$q^\ast_{j+1}$')
zlabel('$q^\ast_{j+2}$')
grid on;
legend;
title('State Space Mistuned system')

% Poincare tuned system
figure(11)
scatter(Q(1,1:sol.N_P:end),Q(2,1:sol.N_P:end),25,...
    'MarkerEdgeColor','k','MarkerFaceColor',color.ies,...
    'DisplayName','$j=0$');
hold on;
box on;
axis tight;
xlabel('$q_j$')
ylabel('$q_{j+1}$')
legend;
title('Poincar\''e Map Tuned system')

% Poincare Mistuned system
figure(12)
scatter(Q_mt(1,1:sol.N_P:end),Q_mt(2,1:sol.N_P:end),25,...
    'MarkerEdgeColor','k','MarkerFaceColor',color.ies,...
    'DisplayName','$j=0$');
hold on;
box on;
axis tight;
xlabel('$q^\ast_j$')
ylabel('$q^\ast_{j+1}$')
legend;
title('Poincar\''e Map Mistuned system')
