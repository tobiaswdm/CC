%% Initialize variables
% Amplitude magnification factors

% Linear case
A_ref = zeros(1, simsetup.ConvergenceMCS.N_MCS);
% Nonlinear case
A = zeros(1,simsetup.ConvergenceMCS.N_MCS);

A_dev_perm = zeros(simsetup.ConvergenceMCS.Shuffle_MCS,...
    simsetup.ConvergenceMCS.N_MCS);
A_ref_dev_perm = zeros(simsetup.ConvergenceMCS.Shuffle_MCS,...
    simsetup.ConvergenceMCS.N_MCS);

% 95% estimate
steps = 48:24:(0.5*simsetup.ConvergenceMCS.N_MCS);
A_95 = zeros(simsetup.ConvergenceMCS.Shuffle_MCS,...
              length(steps));
A_ref_95 = zeros(simsetup.ConvergenceMCS.Shuffle_MCS,...
              length(steps));

%% Build system

% Build nominal system
[sys,exc] = BuildSystem(sys,exc,'tuned');

% Determine performance of tuned system with absorber

% Frequency range for stepping
r_range = [sys.r_k(exc.k+1), sys.r_k_noabs(exc.k+1)];
r_range = simsetup.ConvergenceMCS.r_scale.*r_range;
r_steps = linspace(r_range(1),r_range(2),...
                   simsetup.ConvergenceMCS.N_rSteps);

% Resonance amplitude
[qhat_tuned, qhat_tuned_std, ~] = FindResonance(sys,sol,exc,r_steps,'tuned');

parfor (k = 1:simsetup.ConvergenceMCS.N_MCS,...
                sol.N_Workers)
    disp(['Simulation ' num2str(k) ' of ' num2str(simsetup.ConvergenceMCS.N_MCS)])
    % Build mistuned system
    [sys_mt,exc_mt] = BuildSystem(sys,exc,'mistuned');

    % Linear mistuning analysis

    % Fixed absorbers
    % Frequency response function of max amplitude
    r_lin = linspace(0.99*sys_mt.r_k_mt(1),...
                       1.01*sys_mt.r_k_mt(end),3000);
    q = max(abs(ComputeLinearResponse(r_lin,sys_mt,exc,...
                                'mistuned','fixed_absorbers')),[],1);
    % Extract overall maximum and corresponding frequency
    [q_max,i_max] = max(q,[],2);
    r_range = r_lin(i_max);
    
    % Store linear amplitude magnification 
    A_ref(k) = q_max/sys_mt.qref;

    % Removed absorbers
    % Frequency response function of max amplitude
    r_lin = linspace(0.99*sys_mt.r_k_noabs_mt(1),...
                       1.01*sys_mt.r_k_noabs_mt(end),3000);
    q = max(abs(ComputeLinearResponse(r_lin,sys_mt,exc,...
                                'mistuned','removed_absorbers')),[],1);

    % Extract overall maximum and corresponding frequency
    [~,i_max] = max(q,[],2);
    r_range = horzcat(r_range,r_lin(i_max));

    % Scale frequency range for stepping
    r_range = simsetup.ConvergenceMCS.r_scale.*r_range;
    r_steps = linspace(r_range(1),r_range(2),...
               simsetup.ConvergenceMCS.N_rSteps);

    % Determine resonance amplitude
    [qhat_mt, qhat_mt_std, ~] = ....
    FindResonance(sys_mt,sol,exc,r_steps,'mistuned');
    
    A(k) = qhat_mt/qhat_tuned;

end

% Save Amplitude magnification factors
save([savepath 'A.mat'],'A')
save([savepath 'A_ref.mat'],'A_ref')

%% Mean squared convergence
for i = 1:simsetup.ConvergenceMCS.Shuffle_MCS
    disp(i)
    % Get empirical 99% Intervals with Bootstraping
    parfor j = 1:length(steps)
        
        % Random sample from initial vector with resampling enabled
        ii = randsample(simsetup.ConvergenceMCS.N_MCS,steps(j),true);

        A_95(i,j) = quantile(A(ii),0.95);

        A_ref_95(i,j) = quantile(A_ref(ii),0.95);

    end

end

% Normalize with quantile of full dataset
A_95 = A_95/quantile(A,0.95);
A_ref_95 = A_ref_95/quantile(A_ref,0.95);

A_95_std = std(A_95,0,1);
A_ref_95_std = std(A_ref_95,0,1);

%% Plotting

figure(1)
s=scatterhistogram(A,A_ref,'HistogramDisplayStyle','bar','MarkerSize',1);
s.Color = {color.ies};
xlabel('$A$')
ylabel('$A_\mathrm{ref}$')
title(['Correlation: ' num2str(corr(A',A_ref'))])
savefig([savepath 'A_A_ref_relation.fig'])


tolerance = 5e-2;
[~,ii] = min(abs(3*A_95_std-tolerance));

figure(2)
plot(steps,max(abs(A_95-1),[],1),'--','LineWidth',1.5,'Color',color.show,...
    'DisplayName','MCS Max')
hold on;
plot(steps,3*A_95_std,'LineWidth',1.5,'Color',color.ies,...
    'DisplayName','MCS')
yline(tolerance,'--k','DisplayName','Tolerance')
xline(steps(ii),'--g','DisplayName','Chosen $N_\mathrm{MCS}$')
box on;
axis tight;
xlabel('$N_\mathrm{MCS}$')
ylabel('Error')
title('95\% quantile of $A$')
set(gca,'YScale','log')
legend;
savefig([savepath 'A_95.fig'])

figure(3)
plot(steps,max(abs(A_ref_95-1),[],1),'--','LineWidth',1.5,'Color',color.show,...
    'DisplayName','MCS')
hold on;
plot(steps,3*A_ref_95_std,'LineWidth',1.5,'Color',color.ies,...
    'DisplayName','MCS')
yline(5e-2,'--k','DisplayName','Tolerance')
xline(steps(ii),'--g','DisplayName','Chosen $N_\mathrm{MCS}$')
box on;
axis tight;
xlabel('$N_\mathrm{MCS}$')
ylabel('Error')
title('95\% quantile of $A_\mathrm{ref}$')
set(gca,'YScale','log')
savefig([savepath 'A_ref_95.fig'])

xx = linspace(1.1*min(A_95(:,ii)-1),1.1*max(A_95(:,ii)-1),2000);
pn = fitdist(A_95(:,ii)-1,'Normal');

figure(4)
histogram(A_95(:,ii)-1,'Normalization','pdf','FaceAlpha',0.2)
hold on;
plot(xx,pn.pdf(xx),'LineWidth',1.5)
xline(3*A_95_std(ii),'--k','LineWidth',1.5)
xline(-3*A_95_std(ii),'--k','LineWidth',1.5)
xlabel('Error at chosen $N_\mathrm{MCS}$')
ylabel('PDF')
title('95\% quantile of $A_\mathrm{ref}$')
savefig([savepath 'A_95_hist.fig'])


% CDF
[F,A_plot] = ecdf(A);
[Fmcs,Amcs_plot] = ecdf(A(1:steps(ii)));

figure(5)
plot(Amcs_plot,Fmcs,'-','LineWidth',2,'Color','g','DisplayName','MCS');
hold on;
plot(A_plot,F,'--k','LineWidth',2,'DisplayName','Full');
xlabel('$A$')
ylabel('$F(A)$')
box on;
axis tight;
legend;
title('CDF')
savefig([savepath 'A_cdf.fig'])

% CDF
[F,A_plot] = ecdf(A_ref);
[Fmcs,Amcs_plot] = ecdf(A_ref(1:steps(ii)));

figure(6)
plot(Amcs_plot,Fmcs,'-','LineWidth',2,'Color','g','DisplayName','MCS');
hold on;
plot(A_plot,F,'--k','LineWidth',2,'DisplayName','Full');
xlabel('$A_\mathrm{ref}$')
ylabel('$F\left(A_\mathrm{ref}\right)$')
box on;
axis tight;
legend;
title('CDF')
savefig([savepath 'A_ref_cdf.fig'])