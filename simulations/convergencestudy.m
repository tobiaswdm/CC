%% Build linear system
[K_lin,Phi,omega] = linear_system(sys);
[M,C,K,f,WN] = buildsystem(sys,exc,mist,Phi,omega,'tuned');

% Dispersion Diagram
k = 0:floor(sys.N_s/2);
figure(1);
plot(k,sqrt(1+4*sys.kappa_c*sin(k*pi/sys.N_s).^2),'o','Color',color.reference,'MarkerSize',8,'LineWidth',1.5)
hold on;
plot(linspace(0,floor(sys.N_s/2),1000),sqrt(1+4*sys.kappa_c*sin(linspace(0,floor(sys.N_s/2),1000)*pi/sys.N_s).^2),'LineWidth',1.5,'Color',color.reference,'MarkerSize',8)
box on;
title('Dispersion Diagram')
xlabel('Wave number - $k$')
ylabel('Eigenfrequency ratio - $r_k$')
xlim([0 floor(sys.N_s/2)])
ylim([0.99 1.01*omega(sys.N_s)])
savefig([savepath 'dispersion_diagram.fig'])

%% Reference response

% Limit case fixed absorber
r_fixed = sqrt(1+4*sys.kappa_c*sin(exc.k*pi/sys.N_s).^2);
[eta_fixed,~] = calculatelinearresponse(eye(sys.N_s),C,K,f,r_fixed);
q_fixed = Phi*eta_fixed;
a_fixed = abs(q_fixed);

% Limit case absorber removed
r_removed = sqrt(1/(1-sys.epsilon_a)+4*(sys.kappa_c/(1-sys.epsilon_a))*sin(exc.k*pi/sys.N_s).^2);
[eta_removed,eta_dot_removed] = calculatelinearresponse(M,C,K,f,r_removed);
% Clearance normalized by amplitude of fixed absorber case
% This value is the same in each sector for the linear tuned case
g_ref = a_fixed(1);

%% Set up convergence study

% Conservative estimation of the length of the transient response based on
% the deacy of the homogenous solution of the fundamental mode in the
% case of a fixed absorber that should be decayed by 95%
tau_decay = log(100)/sys.D(1);
Tau_fixed = 2*pi/r_fixed;
Tau_removed = 2*pi/r_removed;

% Number of excitation periods
NT = 50:50:500;

% Time step with respect to excitation frequency
N_Sample = 1000;
NP = 1000 * (2).^(0:10);

dtau_fixed = Tau_fixed./NP;
dtau_removed = Tau_removed./NP;
NP_transient_fixed = NP(1)*ceil(tau_decay/Tau_fixed);
NP_transient_removed = NP(1)*ceil(tau_decay/Tau_removed);

N_Save = NP/N_Sample;

% Initialize matrices
AMEAN_AP = zeros(length(NP),length(NT));
AMAX_AP = zeros(length(NP),length(NT));

AMEAN_SMR = zeros(length(NP),length(NT));
AMAX_SMR = zeros(length(NP),length(NT));

xi = linspace(-a_fixed(1),a_fixed(1),timesim.N_Sample);
dist_ap = zeros(length(NP),timesim.N_Sample);
dist_smr = zeros(length(NP),timesim.N_Sample);

q0 = 1e-4 * a_fixed(1) * (2*rand(sys.N_s,1)-1);
eta0ap = Phi\q0;
eta0smr = Phi\q0;
etadot0ap = zeros(sys.N_s,1);
etadot0smr = zeros(sys.N_s,1);
g_scale_ap = 0.1;
qa0ap = -0.99*g_scale_ap*g_ref + q0;
qadot0ap = zeros(sys.N_s,1);
g_scale_smr = 0.4;
qa0smr = -0.99*g_scale_smr*g_ref + q0;
qadot0smr = zeros(sys.N_s,1);


for i = 1:length(NP)
    
    disp(['Sample rate ' num2str(i) ' of ' num2str(length(NP))])
    systemp = sys;
    %% Almost periodic regime
    systemp.g_scale = g_scale_ap;
    [Ma,g,WNA] = buildabsorber(g_ref,systemp,'tuned');
    
    excap = exc;
    excap.harmonic.r = r_fixed;
    % Time simulation
    if i==1
        [ETAAP,QAAP,ETADOTAP,UAAP,TAUAP] = time_moreau(M,Ma,C,K,eta0ap,qa0ap,etadot0ap,qadot0ap,excap,f,WN,WNA,g,systemp.eN,sol,dtau_fixed(i),...
            NP_transient_fixed,NT(end)*NP(i),N_Save(i));
    else
        [ETAAP,QAAP,ETADOTAP,UAAP,TAUAP] = time_moreau(M,Ma,C,K,eta0ap,qa0ap,etadot0ap,qadot0ap,excap,f,WN,WNA,g,systemp.eN,sol,dtau_fixed(i),...
            0,NT(end)*NP(i),N_Save(i));
    end
    QAP = Phi*ETAAP;
    
    eta0ap = ETAAP(:,end);
    etadot0ap = ETADOTAP(:,end);
    qa0ap = QAAP(:,end);
    qadot0ap = UAAP(:,end);
    
    %% SMR regime
    systemp.g_scale = g_scale_smr;
    [Ma,g,WNA] = buildabsorber(g_ref,systemp,'tuned');
    
    excsmr = exc;
    excsmr.harmonic.r = r_removed;
    % Time simulation
    if i==1
        [ETASMR,QASMR,ETADOTSMR,UASMR,TAUSMR] = time_moreau(M,Ma,C,K,eta0smr,qa0smr,etadot0smr,qadot0smr,excsmr,f,WN,WNA,g,systemp.eN,sol,dtau_removed(i),...
        NP_transient_removed,NT(end)*NP(i),N_Save(i));

    else
        [ETASMR,QASMR,ETADOTSMR,UASMR,TAUSMR] = time_moreau(M,Ma,C,K,eta0smr,qa0smr,etadot0smr,qadot0smr,excsmr,f,WN,WNA,g,systemp.eN,sol,dtau_removed(i),...
            0,NT(end)*NP(i),N_Save(i));
    end
    QSMR = Phi*ETASMR;
    
    dist_ap(i,:) = pdf1(QAP(1,:),xi);
    dist_smr(i,:) = pdf1(QSMR(1,:),xi);
    
    eta0smr = ETASMR(:,end);
    etadot0smr = ETADOTSMR(:,end);
    qa0smr = QASMR(:,end);
    qadot0smr = UASMR(:,end);
    
    MEANAP = zeros(systemp.N_s,length(NT));
    MEANSMR = zeros(systemp.N_s,length(NT));
    MAXAP = zeros(systemp.N_s,length(NT));
    MAXSMR = zeros(systemp.N_s,length(NT));
    
    for j = 1:length(NT)
        [MEANAP(:,j),~,MAXAP(:,j)] = meanmaxamp(QAP(:,2:(N_Sample*NT(j)+1)),N_Sample);
        [MEANSMR(:,j),~,MAXSMR(:,j)] = meanmaxamp(QSMR(:,2:(N_Sample*NT(j)+1)),N_Sample);
    end
    
    AMEAN_AP(i,:) = mean(MEANAP,1);
    AMAX_AP(i,:) = max(MAXAP,[],1);
    
    AMEAN_SMR(i,:) = mean(MEANSMR,1);
    AMAX_SMR(i,:) = max(MAXSMR,[],1);
    disp(['Finished sample rate ' num2str(i) ' of ' num2str(length(NP))])
    
    if i == length(1)
        figure(2);
        plot(TAUAP,QAP(1,:))
        hold on;
        xlim([TAUAP(1) TAUAP(end)])
        yline(AMEAN_AP(i,end),'--r','LineWidth',1.5)
        yline(AMAX_AP(i,end),'-r','LineWidth',1.5)
        xlabel('$\tau$')
        ylabel('$q_j$')
        %savefig([savepath 'time_ap.fig'])
        tendap = TAUAP(end);
        
        figure(3);
        plot(TAUSMR,QSMR(1,:))
        hold on;
        xlim([TAUSMR(1) TAUSMR(end)])
        yline(AMEAN_SMR(i,end),'--r','LineWidth',1.5)
        yline(AMAX_SMR(i,end),'-r','LineWidth',1.5)
        xlabel('$\tau$')
        ylabel('$q_j$')
        %savefig([savepath 'time_smr.fig'])
        tendsmr=TAUSMR(end);
    else
        figure(2);
        plot(TAUAP+tendap,QAP(1,:))
        hold on;
        xlim([-inf inf])
        %yline(AMEAN_AP(i,end),'--r','LineWidth',1.5)
        %yline(AMAX_AP(i,end),'-r','LineWidth',1.5)
        xlabel('$\tau$')
        ylabel('$q_j$')
        %savefig([savepath 'time_ap.fig'])
        tendap = TAUAP(end)+tendap;
        
        figure(3);
        plot(TAUSMR+tendsmr,QSMR(1,:))
        hold on;
        xlim([-inf inf])
        %yline(AMEAN_SMR(i,end),'--r','LineWidth',1.5)
        %yline(AMAX_SMR(i,end),'-r','LineWidth',1.5)
        xlabel('$\tau$')
        ylabel('$q_j$')
        %savefig([savepath 'time_smr.fig'])
        tendsmr=TAUSMR(end)+tendsmr;
    end
end

figure;
subplot(2,1,1)
hold on;
plot(NP(1:end-1),abs(AMEAN_AP(1:end-1,end)/AMEAN_AP(end,end)-1),'-+','LineWidth',1.5,'Color',color.ies)
title('Mean Amplitude - Almost periodic')
xlabel('Points per excitation period')
ylabel('Deviation from highest sample rate')
box on;
xlim([NP(1), NP(end-1)])
%ylim([0 0.05])
set(gca,'YScale','lin')
set(gca,'XScale','log')
subplot(2,1,2)
hold on;
plot(NT(1:end-1),abs(AMEAN_AP(end,1:end-1)/AMEAN_AP(end,end)-1),'-+','LineWidth',1.5,'Color',color.ies)
xlabel('Number of excitation periods')
ylabel('Deviation from longest simulation')
box on;
xlim([NT(1), NT(end-1)])
%ylim([0 0.05])
set(gca,'YScale','lin')
set(gca,'XScale','lin')
savefig([savepath 'meanamp_ap.fig'])


figure;
subplot(2,1,1)
hold on;
plot(NP(1:end-1),abs(AMAX_AP(1:end-1,end)'/AMAX_AP(end,end)-1),'-+','LineWidth',1.5,'Color',color.ies)
set(gca,'YScale','lin')
set(gca,'XScale','log')
title('Maximum Amplitude - Almost periodic')
xlabel('Points per excitation period')
ylabel('Deviation from highest sample rate')
box on;
xlim([NP(1), NP(end-1)])
%ylim([0 0.05])
subplot(2,1,2)
hold on;
loglog(NT(1:end-1),abs(AMAX_AP(end,1:end-1)/AMAX_AP(end,end)-1),'-+','LineWidth',1.5,'Color',color.ies)
xlabel('Number of excitation periods')
ylabel('Deviation from longest simulation')
box on;
xlim([NT(1), NT(end-1)])
%ylim([0 0.05])
set(gca,'YScale','lin')
set(gca,'XScale','lin')
savefig([savepath 'maxamp_ap.fig'])


figure;
subplot(2,1,1)
hold on;
plot(NP(1:end-1),abs(AMEAN_SMR(1:end-1,end)'/AMEAN_SMR(end,end)-1),'-+','LineWidth',1.5,'Color',color.ies)
title('Mean Amplitude - SMR')
xlabel('Points per excitation period')
ylabel('Deviation from highest sample rate')
box on;
xlim([NP(1), NP(end-1)])
%ylim([0 0.05])
set(gca,'YScale','lin')
set(gca,'XScale','log')
subplot(2,1,2)
hold on;
plot(NT(1:end-1),abs(AMEAN_SMR(end,1:end-1)/AMEAN_SMR(end,end)-1),'-+','LineWidth',1.5,'Color',color.ies)
xlabel('Number of excitation periods')
ylabel('Deviation from longest simulation')
box on;
xlim([NT(1), NT(end-1)])
%ylim([0 0.05])
set(gca,'YScale','lin')
set(gca,'XScale','lin')
savefig([savepath 'meanamp_smr.fig'])

figure;
subplot(2,1,1)
hold on;
plot(NP(1:end-1),abs(AMAX_SMR(1:end-1,end)'/AMAX_SMR(end,end)-1),'-+','LineWidth',1.5,'Color',color.ies)
title('Maximum Amplitude - SMR')
xlabel('Points per excitation period')
ylabel('Deviation from highest sample rate')
box on;
xlim([NP(1), NP(end-1)])
%ylim([0 0.05])
set(gca,'YScale','lin')
set(gca,'XScale','log')
subplot(2,1,2)
hold on;
plot(NT(1:end-1),abs(AMAX_SMR(end,1:end-1)/AMAX_SMR(end,end)-1),'-+','LineWidth',1.5,'Color',color.ies)
xlabel('Number of excitation periods')
ylabel('Deviation from longest simulation')
box on;
xlim([NT(1), NT(end-1)])
set(gca,'YScale','lin')
set(gca,'XScale','lin')
savefig([savepath 'maxamp_smr.fig'])

figure;
plot(xi,dist_ap)
hold on;
box on;
xlim([xi(1) xi(end)])
legend(string(NP))
xlabel('$q_0$')
ylabel('PDF')
title('Almost periodic regime')
savefig([savepath 'pdf_ap.fig'])
cols = hot(8);
figure;
hold on;
for i = 1:6
    plot(xi/a_fixed(1),dist_smr(i,:),'Color',cols(7-(i-1),:),'LineWidth',1.5)
end
box on;
xlim([-0.15 0.15])
legend(string(NP(1:6)))
xlabel('$q_0$')
ylabel('PDF')
title('SMR regime')
savefig([savepath 'pdf_smr.fig'])

save([savepath 'AMEAN_SMR.mat'],'AMEAN_SMR')
save([savepath 'AMAX_SMR.mat'],'AMAX_SMR')
save([savepath 'AMEAN_AP.mat'],'AMEAN_AP')
save([savepath 'AMAX_AP.mat'],'AMAX_AP')
save([savepath 'NP.mat'],'NP')
save([savepath 'dist_ap.mat'],'dist_ap')
save([savepath 'dist_smr.mat'],'dist_smr')
save([savepath 'xi.mat'],'xi')