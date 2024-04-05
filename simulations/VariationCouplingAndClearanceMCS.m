%% Initialize variables
% Amplitude magnification factors

% Linear case
A_ref = zeros(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c, ...
              simsetup.VariationCouplingAndClearanceMCS.N_MCS);
% Nonlinear case
A = zeros(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c, ...
          simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale,...
          simsetup.VariationCouplingAndClearanceMCS.N_MCS);
A_max = zeros(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c, ...
              simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale);
delta_omega_A = zeros(sys.N_s,simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c, ...
                 simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale,...
                 simsetup.VariationCouplingAndClearanceMCS.N_MCS);
delta_omega_A_max = zeros(sys.N_s,simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c, ...
                 simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale);
% 95% quantile Estimate 
A_95 = zeros(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c, ...
             simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale);
A_ref_95 = zeros(1,simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c);

% Actual amplitudes
qhat_tuned = zeros(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c, ...
                simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale);
qhat_tuned_std = zeros(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c, ...
                simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale);
% Number of significant impacts per excitation period 
N_sipp_tuned = zeros(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c, ...
                simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale);

qhat_mt = zeros(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c, ...
                simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale, ...
                simsetup.VariationCouplingAndClearanceMCS.N_MCS);
qhat_std_mt = zeros(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c, ...
                    simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale, ...
                    simsetup.VariationCouplingAndClearanceMCS.N_MCS);
% Number of significant impacts per excitation period 
N_sipp_mt = zeros(sys.N_s,simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c, ...
                  simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale, ...
                  simsetup.VariationCouplingAndClearanceMCS.N_MCS);

% Response types
Resp_type_tuned = nan(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c, ...
                  simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale);
Resp_type_mt = nan(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c, ...
                  simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale, ...
                  simsetup.VariationCouplingAndClearanceMCS.N_MCS);

% Maximum linear Amplitude magnification factor
A_ref_max = zeros(1,simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c);
delta_omega_Aref = zeros(sys.N_s,simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c,...
                         simsetup.VariationCouplingAndClearanceMCS.N_MCS);
delta_omega_Aref_max = zeros(sys.N_s,simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c);


% Determine couplings strengts
switch simsetup.VariationCouplingAndClearanceMCS.Scaling_kappa_c
    case 'logarithmic'
        kappa_c = logspace(log10(simsetup.VariationCouplingAndClearanceMCS.Range_kappa_c(1)),...
                           log10(simsetup.VariationCouplingAndClearanceMCS.Range_kappa_c(2)),...
                           simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c);
    case 'linear'
        kappa_c = linspace(simsetup.VariationCouplingAndClearanceMCS.Range_kappa_c(1),...
                           simsetup.VariationCouplingAndClearanceMCS.Range_kappa_c(2),...
                           simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c);
    otherwise
        error('Case not defined.')
end

% Determine normalized clearances
switch simsetup.VariationCouplingAndClearanceMCS.Scaling_GammaScale
    case 'logarithmic'
        Gamma_Scale = logspace(log10(simsetup.VariationCouplingAndClearanceMCS.Range_GammaScale(1)),...
                           log10(simsetup.VariationCouplingAndClearanceMCS.Range_GammaScale(2)),...
                           simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale);
    case 'linear'
        Gamma_Scale = linspace(simsetup.VariationCouplingAndClearanceMCS.Range_GammaScale(1),...
                           simsetup.VariationCouplingAndClearanceMCS.Range_GammaScale(2),...
                           simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale);
    otherwise
        error('Case not defined.')
end

save([savepath 'Gamma_Scale.mat'],'Gamma_Scale')
save([savepath 'kappa_c.mat'],'kappa_c')



%% Vary nominal parameters

for i = 1:simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c
    disp(['Coupling strength ' num2str(i) ' of ' ...
        num2str(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c)])

    % Set coupling strength
    sys.kappa_c = kappa_c(i);

    for j = 1:simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale
        disp(['Clearance ' num2str(j) ' of ' ...
        num2str(simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale)])

    % Set nominal clearance
    sys.Gamma_Scale = Gamma_Scale(j);

    % Build nominal system
    [sys,exc] = BuildSystem(sys,exc,'tuned');

    % Determine performance of tuned system with absorber

    % Frequency range for stepping
    r_range = [sys.r_k(exc.k+1), sys.r_k_noabs(exc.k+1)];
    r_range = simsetup.VariationCouplingAndClearanceMCS.r_scale.*r_range;
    r_steps = linspace(r_range(1),r_range(2),...
                       simsetup.VariationCouplingAndClearanceMCS.N_rSteps);
    
    % Resonance amplitude
    [qhat_tuned(i,j), qhat_tuned_std(i,j), N_sipp_tuned(i,j)] = ....
    FindResonance(sys,sol,exc,r_steps,'tuned');
    
    % Response type
    if N_sipp_tuned(i,j)>=1.98
        Resp_type_tuned(i,j) = 0; % GSAPR
    else
        Resp_type_tuned(i,j) = 2; % SMR
    end

    
    disp('Starting MCS loop')
        parfor (k = 1:simsetup.VariationCouplingAndClearanceMCS.N_MCS,...
                sol.N_Workers)

            
            % Build mistuned system
            [sys_mt,exc_mt] = BuildSystem(sys,exc,'mistuned');

            % Save mistuning configuration
            delta_omega_A(:,i,j,k) = sys_mt.delta_omega;
            
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
            
            % Store linear amplitude magnification only for first
            % nominal clearance (only depends on coupling in linear case)
            if j == 1
                A_ref(i,k) = q_max/sys_mt.qref;
                delta_omega_Aref(:,i,k) = sys_mt.delta_omega;
            end
            
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
            r_range = simsetup.VariationCouplingAndClearanceMCS.r_scale.*r_range;
            r_steps = linspace(r_range(1),r_range(2),...
                       simsetup.VariationCouplingAndClearanceMCS.N_rSteps);

            % Determine resonance amplitude
            [qhat_mt_temp, qhat_mt_std_temp, N_sipp_mt_temp] = ....
            FindResonance(sys_mt,sol,exc,r_steps,'mistuned');
            qhat_mt(i,j,k) = qhat_mt_temp;
            qhat_mt_std(i,j,k) = qhat_mt_std_temp;
            N_sipp_mt(:,i,j,k) = N_sipp_mt_temp;
            
            % Save backup
            parsave(savepath_backup,i,j,k,qhat_tuned, qhat_tuned_std,...
                N_sipp_tuned, qhat_mt_temp, qhat_mt_std_temp, N_sipp_mt_temp)

            % Compute amplitude magnification factor
            A(i,j,k) = qhat_mt_temp/qhat_tuned(i,j);

            % Response type
            Resp_type_mt(i,j,k) = ResponseType(N_sipp_mt_temp);
        end

        

        % Maximum amplification
        [A_max(i,j), k_max] = max(A(i,j,:));
        delta_omega_A_max(:,i,j) = delta_omega_A(:,i,j,k_max);

        % 95% Interval
        A_95(i,j) = quantile(squeeze(A(i,j,:)),0.95);

        % Save data
        save([savepath 'A.mat'],'A')
        save([savepath 'A_max.mat'],'A_max')
        save([savepath 'A_95.mat'],'A_95')
        save([savepath 'A_ref.mat'],'A_ref')
        save([savepath 'delta_omega_Aref_max.mat'],'delta_omega_Aref_max')
        save([savepath 'delta_omega_A_max.mat'],'delta_omega_A_max')
        save([savepath 'Resp_type_mt.mat'],'Resp_type_mt')
        save([savepath 'Resp_type_tuned.mat'],'Resp_type_tuned')
        save([savepath 'qhat_mt.mat'],'qhat_mt')
        save([savepath 'qhat_mt_std.mat'],'qhat_mt_std')

    end

    % Find system with largest linear Amplitude magnification factor
    [A_ref_max(i), k_max] = max(A_ref(i,:),[],2);
    delta_omega_Aref_max(:,i) = delta_omega_Aref(:,i,k_max);
    save([savepath 'A_ref_max.mat'],'A_ref_max')

    % 95% Interval
    A_ref_95(i) = quantile(A_ref(i,:),0.95);
    save([savepath 'A_ref_95.mat'],'A_ref_95')

end

% Save data
save([savepath 'A.mat'],'A')
save([savepath 'A_max.mat'],'A_max')
save([savepath 'A_95.mat'],'A_95')
save([savepath 'A_ref.mat'],'A_ref')
save([savepath 'A_ref_max.mat'],'A_ref_max')
save([savepath 'A_ref_95.mat'],'A_ref_95')
save([savepath 'delta_omega_Aref_max.mat'],'delta_omega_Aref_max')
save([savepath 'delta_omega_A_max.mat'],'delta_omega_A_max')
save([savepath 'Resp_type_mt.mat'],'Resp_type_mt')
save([savepath 'Resp_type_tuned.mat'],'Resp_type_tuned')
save([savepath 'qhat_mt.mat'],'qhat_mt')
save([savepath 'qhat_mt_std.mat'],'qhat_mt_std')


figure(1)
s=scatterhistogram(squeeze(A(1,1,:)),A_ref(1,:),'HistogramDisplayStyle',...
    'bar','MarkerSize',2);
s.Color = {color.ies};
xlabel('$A$')
ylabel('$A_\mathrm{ref}$')
title(['Correlation: ' num2str(corr(A_ref(1,:)',squeeze(A(1,1,:))))])
savefig([savepath 'A_A_ref_relation.fig'])


if (simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale>1 && ...
    simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c>1)

    [XX,YY] = meshgrid(Gamma_Scale,kappa_c);
    % Max A
    figure(2);
    surf(XX,YY,A_max,'EdgeAlpha',0)
    hold on;
    title('Maximum from MCSs')
    box on;
    colormap turbo
    xlabel('$\Gamma/\hat{q}_\mathrm{ref}$')
    ylabel('$\kappa_\mathrm{c}$')
    zlabel('$\mathrm{max}\left\{A \right\}$')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    xlim([simsetup.VariationCouplingAndClearanceMCS.Range_GammaScale(1),...
         simsetup.VariationCouplingAndClearanceMCS.Range_GammaScale(2)])
    ylim([simsetup.VariationCouplingAndClearanceMCS.Range_kappa_c(1),...
         simsetup.VariationCouplingAndClearanceMCS.Range_kappa_c(2)])
    zlim([-inf inf])
    savefig([savepath 'A_max.fig'])
    
    % 95% Interval
    figure(3);
    surf(XX,YY,A_95,'EdgeAlpha',0)
    hold on;
    title('95\% percentile w/ absorbers')
    box on;
    colormap turbo
    xlabel('$\Gamma/\hat{q}_\mathrm{ref}$')
    ylabel('$\kappa_\mathrm{c}$')
    zlabel('$A | \mathcal{F} (A) = 95\, \%$')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    xlim([simsetup.VariationCouplingAndClearanceMCS.Range_GammaScale(1),...
         simsetup.VariationCouplingAndClearanceMCS.Range_GammaScale(2)])
    ylim([simsetup.VariationCouplingAndClearanceMCS.Range_kappa_c(1),...
         simsetup.VariationCouplingAndClearanceMCS.Range_kappa_c(2)])
    zlim([-inf inf])
    savefig([savepath 'A_95.fig'])
    
    
    figure(4);
    plot(kappa_c,A_ref_95,'LineWidth',1.5,'Color',color.reference)
    hold on;
    title('95\% percentile w/o absorbers')
    xlabel('$\kappa_\mathrm{c}$')
    ylabel('$A_\mathrm{ref} | \mathcal{F} (A_\mathrm{ref}) = 95\, \%$')
    xlim([simsetup.VariationCouplingAndClearanceMCS.Range_kappa_c(1),...
         simsetup.VariationCouplingAndClearanceMCS.Range_kappa_c(2)])
    ylim([-inf inf])
    savefig([savepath 'A_ref_95.fig'])
    
    % Max weibull
    figure(5);
    plot(kappa_c,A_ref_max,'LineWidth',1.5,'Color',color.reference)
    hold on;
    title('Maximum w/o absorbers')
    xlabel('$\kappa_\mathrm{c}$')
    ylabel('$\mathrm{max}\left\{A_\mathrm{ref} \right\}$')
    xlim([simsetup.VariationCouplingAndClearanceMCS.Range_kappa_c(1),...
         simsetup.VariationCouplingAndClearanceMCS.Range_kappa_c(2)])
    ylim([-inf inf])
    savefig([savepath 'A_ref_max.fig'])

elseif (simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale==1 && ...
        simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c>1)
    % A max
    figure(2);
    plot(kappa_c,A_max,'LineWidth',1.5,'Color',color.ies)
    hold on;
    title('Maximum from MCSs')
    box on;
    xlabel('$\kappa_\mathrm{c}$')
    ylabel('$\mathrm{max}\left\{A \right\}$')
    set(gca,'XScale','log')
    xlim([simsetup.VariationCouplingAndClearanceMCS.Range_kappa_c(1),...
         simsetup.VariationCouplingAndClearanceMCS.Range_kappa_c(2)])
    ylim([-inf inf])
    savefig([savepath 'A_max.fig'])
    
    % 95% Interval
    figure(3);
    plot(kappa_c,A_95,'LineWidth',1.5,'Color',color.ies)
    hold on;
    title('95\% percentile w/ absorbers')
    box on;
    xlabel('$\kappa_\mathrm{c}$')
    ylabel('$A | \mathcal{F} (A) = 95\, \%$')
    set(gca,'XScale','log')
    xlim([simsetup.VariationCouplingAndClearanceMCS.Range_kappa_c(1),...
         simsetup.VariationCouplingAndClearanceMCS.Range_kappa_c(2)])
    ylim([-inf inf])
    savefig([savepath 'A_95.fig'])
   
    % A ref Weibull
    figure(4);
    plot(kappa_c,A_ref_95,'LineWidth',1.5,'Color',color.reference)
    hold on;
    title('95\% percentile w/o absorbers')
    xlabel('$\kappa_\mathrm{c}$')
    ylabel('$A_\mathrm{ref} | \mathcal{F} (A_\mathrm{ref}) = 95\, \%$')
    xlim([simsetup.VariationCouplingAndClearanceMCS.Range_kappa_c(1),...
         simsetup.VariationCouplingAndClearanceMCS.Range_kappa_c(2)])
    ylim([-inf inf])
    savefig([savepath 'A_ref_95.fig'])

    % A ref max
    figure(5);
    plot(kappa_c,A_ref_max,'LineWidth',1.5,'Color',color.reference)
    hold on;
    title('95\% percentile w/o absorbers')
    xlabel('$\kappa_\mathrm{c}$')
    ylabel('$\mathrm{max}\left\{A_\mathrm{ref} \right\}$')
    xlim([simsetup.VariationCouplingAndClearanceMCS.Range_kappa_c(1),...
         simsetup.VariationCouplingAndClearanceMCS.Range_kappa_c(2)])
    ylim([-inf inf])
    savefig([savepath 'A_ref_max.fig'])

elseif (simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale>1 && ...
        simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c==1)

    % A max
    figure(2);
    plot(Gamma_Scale,A_max,'LineWidth',1.5,'Color',color.ies)
    hold on;
    title('Maximum from MCSs')
    box on;
    xlabel('$\Gamma/\hat{q}_\mathrm{ref}$')
    ylabel('$\mathrm{max}\left\{A \right\}$')
    set(gca,'XScale','log')
    axis tight;
    xlim([simsetup.VariationCouplingAndClearanceMCS.Range_GammaScale(1),...
         simsetup.VariationCouplingAndClearanceMCS.Range_GammaScale(2)])
    savefig([savepath 'A_max.fig'])
    
    % 95% Interval
    figure(3);
    plot(Gamma_Scale,A_95,'LineWidth',1.5,'Color',color.ies)
    hold on;
    title('95\% percentile w/ absorbers')
    box on;
    xlabel('$\Gamma/\hat{q}_\mathrm{ref}$')
    ylabel('$A | \mathcal{F} (A) = 95\, \%$')
    set(gca,'XScale','log')
    xlim([simsetup.VariationCouplingAndClearanceMCS.Range_GammaScale(1),...
         simsetup.VariationCouplingAndClearanceMCS.Range_GammaScale(2)])
    ylim([-inf inf])
    savefig([savepath 'A_95.fig'])
   
    disp(['95% Interval A_ref estimation: ' num2str(A_ref_95)])
    disp(['Maximum A_ref: ' num2str(A_ref_max)])

else
    disp(['95% Interval A estimation: ' num2str(A_95)])
    disp(['Maximum A: ' num2str(A_max)])
    disp(['95% Interval A_ref estimation: ' num2str(A_ref_95)])
    disp(['Maximum A_ref: ' num2str(A_ref_max)])
end


figure(6)
histogram(squeeze(Resp_type_mt(1,1,:)),'Normalization','pdf',...
    'FaceColor',color.ies);
hold on;
set(gca,'XTick',[0 1 2],'XTickLabel',{'GSAPR';'LSAPR';'SMR'})
ylabel('PDF')
title('Solution types')
axis tight;
savefig([savepath 'Resp_type_mt.fig'])

figure(7)
s=scatterhistogram(squeeze(qhat_mt(1,1,:))/sys.qref,A_ref(1,:),'HistogramDisplayStyle','bar');
s.Color = {color.ies};
xlabel('$\mathrm{max}_j \left\{ \hat{q}_j^\ast \right\} / \hat{q}_\mathrm{ref}$')
ylabel('$A_\mathrm{ref}$')
title(['Correlation: ' num2str(corr(A_ref(1,:)',squeeze(qhat_mt(1,1,:))/sys.qref))])
savefig([savepath 'qhat_A_ref_relation.fig'])



function parsave(savepath_backup,i,j,k,qhat_tuned, qhat_tuned_std, N_sipp_tuned, qhat_mt, qhat_std_mt, N_sipp_mt)
save([savepath_backup 'qhat_tuned_' num2str(i) '_' num2str(j) '_' num2str(k) '.mat'],'qhat_tuned')
save([savepath_backup 'qhat_tuned_std_' num2str(i) '_' num2str(j) '_' num2str(k) '.mat'],'qhat_tuned_std')
save([savepath_backup 'N_sipp_tuned_' num2str(i) '_' num2str(j) '_' num2str(k) '.mat'],'N_sipp_tuned')
save([savepath_backup 'qhat_mt_' num2str(i) '_' num2str(j) '_' num2str(k) '.mat'],'qhat_mt')
save([savepath_backup 'qhat_std_mt_' num2str(i) '_' num2str(j) '_' num2str(k) '.mat'],'qhat_std_mt')
save([savepath_backup 'N_sipp_mt_' num2str(i) '_' num2str(j) '_' num2str(k) '.mat'],'N_sipp_mt')
end