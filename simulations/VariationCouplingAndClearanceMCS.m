%% Initialize variables
% Amplitude magnification factors

if sys.sigma_omega ~= 0
    % Linear case if linear mistuning enabled
    A_ref = zeros(...
        simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c, ...
              simsetup.VariationCouplingAndClearanceMCS.N_MCS);

    % 95% quantile Estimate 
    A_ref_95 = zeros(1,simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c);

    % Maximum linear Amplitude magnification factor
    A_ref_max = zeros(...
        1,simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c);
    delta_omega_Aref = zeros(...
        sys.N_s,simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c,...
                     simsetup.VariationCouplingAndClearanceMCS.N_MCS);
    delta_omega_Aref_max = zeros(...
        sys.N_s,simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c);

end

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
A_95 = zeros(...
    simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c, ...
         simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale);


% Actual amplitudes
qref = zeros(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c,1);
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
        Gamma_Scale = logspace(...
            log10(simsetup.VariationCouplingAndClearanceMCS.Range_GammaScale(1)),...
               log10(simsetup.VariationCouplingAndClearanceMCS.Range_GammaScale(2)),...
               simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale);
    case 'linear'
        Gamma_Scale = linspace(...
            simsetup.VariationCouplingAndClearanceMCS.Range_GammaScale(1),...
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
    qref(i) = sys.qref;

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
            if (j == 1) && (sys_mt.sigma_omega~=0)
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
            r_range = simsetup.VariationCouplingAndClearanceMCS.r_scale.*...
                r_range;
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
        save([savepath 'delta_omega_A_max.mat'],'delta_omega_A_max')
        save([savepath 'Resp_type_mt.mat'],'Resp_type_mt')
        save([savepath 'Resp_type_tuned.mat'],'Resp_type_tuned')
        save([savepath 'qhat_mt.mat'],'qhat_mt')
        save([savepath 'qhat_mt_std.mat'],'qhat_mt_std')
        if sys.sigma_omega ~= 0
            save([savepath 'A_ref.mat'],'A_ref')
            save([savepath 'delta_omega_Aref_max.mat'],...
                'delta_omega_Aref_max')
        end

    end
    
    if sys.sigma_omega ~= 0
        % Find system with largest linear Amplitude magnification factor
        [A_ref_max(i), k_max] = max(A_ref(i,:),[],2);
        delta_omega_Aref_max(:,i) = delta_omega_Aref(:,i,k_max);
        save([savepath 'A_ref_max.mat'],'A_ref_max')
    
        % 95% Interval
        A_ref_95(i) = quantile(A_ref(i,:),0.95);
        save([savepath 'A_ref_95.mat'],'A_ref_95')
    end

end

% Save data
save([savepath 'A.mat'],'A')
save([savepath 'A_max.mat'],'A_max')
save([savepath 'A_95.mat'],'A_95')
save([savepath 'delta_omega_A_max.mat'],'delta_omega_A_max')
save([savepath 'Resp_type_mt.mat'],'Resp_type_mt')
save([savepath 'Resp_type_tuned.mat'],'Resp_type_tuned')
save([savepath 'qhat_mt.mat'],'qhat_mt')
save([savepath 'qhat_mt_std.mat'],'qhat_mt_std')
save([savepath 'qhat_tuned.mat'],'qhat_tuned')
save([savepath 'qref.mat'],'qref')

% Vectors containing the indices
i_ind = repelem(1:simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c, ...
       simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale);
j_ind = repmat(1:simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale, ...
       [1, simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c]);

% Plot Scatter Histograms
if sys.sigma_omega ~= 0
    save([savepath 'A_ref.mat'],'A_ref')
    save([savepath 'A_ref_max.mat'],'A_ref_max')
    save([savepath 'A_ref_95.mat'],'A_ref_95')
    save([savepath 'delta_omega_Aref_max.mat'],'delta_omega_Aref_max')
    
    for i = 1:(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c*...
             simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale)
        
        figure(1)
        subplot(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c,...
            simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale,i)
        s=scatterhistogram(squeeze(A(i_ind(i),j_ind(i),:)),A_ref(i_ind(i),:),...
        'HistogramDisplayStyle','bar','MarkerSize',1,'LineStyle','none');
        s.Color = {color.ies};
        s.LineWidth = 0.1;
        xlabel('$A$')
        ylabel('$A_\mathrm{ref}$')
        title(['\rho = ' num2str(round(1000*corr(A_ref(i_ind(i),:)', ...
        squeeze(A(i_ind(i),j_ind(i),:))))/1000)])

        figure(2)
        subplot(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c,...
            simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale,i)
        s=scatterhistogram(squeeze(qhat_mt(i_ind(i),j_ind(i),:))/qref(i_ind(i)) ...
            ,A_ref(i_ind(i),:),...
        'HistogramDisplayStyle','bar','MarkerSize',1,'LineStyle','none');
        s.Color = {color.ies};
        s.LineWidth = 0.1;
        xlabel('$\hat{q}^\ast / \hat{q}_\mathrm{ref}$')
        ylabel('$A_\mathrm{ref}$')
        title(['\rho = ' num2str(round(1000*corr(A_ref(i_ind(i),:)', ...
        squeeze(qhat_mt(i_ind(i),j_ind(i),:))))/1000)])
    end

end


for i = 1:(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c*...
         simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale)
    
    [E,x] = ecdf(squeeze(A(i_ind(i),j_ind(i),:)));

    figure(3)
    subplot(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c,...
        simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale,i)
    hold on;
    yyaxis right;
    stairs(x,E,'--k','LineWidth',1.5)
    ylabel('CDF')
    yyaxis left;
    histogram(squeeze(A(i_ind(i),j_ind(i),:)),...
    'Normalization','pdf','FaceColor',color.ies,'LineStyle','none');
    box on;
    xlabel('$A$')
    ylabel('PDF')
    axis tight
    
    [E,x] = ecdf(squeeze(qhat_mt(i_ind(i),j_ind(i),:))/qref(i_ind(i)));

    figure(4)
    subplot(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c,...
        simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale,i)
    hold on;
    yyaxis right;
    stairs(x,E,'--k','LineWidth',1.5)
    ylabel('CDF')
    yyaxis left;
    histogram(squeeze(qhat_mt(i_ind(i),j_ind(i),:))/qref(i_ind(i)),...
    'Normalization','pdf','FaceColor',color.ies,'LineStyle','none');
    box on;
    xlabel('$\hat{q}^\ast / \hat{q}_\mathrm{ref}$')
    ylabel('PDF')
    axis tight

end


figure(3)
savefig([savepath 'A_PDF.fig'])

figure(4)
savefig([savepath 'qhat_PDF.fig'])


% Plot response types
figure(5)
for i = 1:(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c*...
             simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale)

    subplot(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c,...
            simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale,i)
    histogram(squeeze(Resp_type_mt(i_ind(i),j_ind(i),:)),'Normalization', ...
        'pdf','FaceColor',color.ies,'FaceAlpha',1,'LineWidth',1);
    hold on;
    set(gca,'XTick',[0 1 2],'XTickLabel',{'GSR';'LSR';'SMR'})
    ylabel('PDF')
    xlim([-0.5 2.5])
    ylim([0 1])
    box on;      
end
savefig([savepath 'Resp_type_mt.fig'])


% Colorscale
if simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale==3
    cScale = [myColors('orange');myColors('cyan');myColors('black')];
else
    cScale = 1-hot(round( ...
    1.8*simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale));
    cScale = flipud( ...
        cScale(1:simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale,:));
end

figure(6)
hold on;
for i = 1:simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale
    yyaxis left;
    plot(kappa_c,A_95(:,i),'-o','LineWidth',1.5,'Color',cScale(i,:), ...
        'MarkerFaceColor',cScale(i,:))
    box on;
    yyaxis right;
    plot(kappa_c,(A_95(:,i).*qhat_tuned(:,i))./qref ...
        ,'--+','LineWidth',1.5,'Color',cScale(i,:))
end
set(gca,'XScale','log')
xlabel('$\kappa_\mathrm{c}$')
yyaxis left
ylabel('$A_{95}$')
axis tight;
yyaxis right
ylabel('$(\hat{q}^\ast / \hat{q}_\mathrm{ref})_{95}$')
savefig([savepath 'A_95.fig'])

figure(7)
hold on;
for i = 1:simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale
    yyaxis left;
    plot(kappa_c,A_max(:,i),'-o','LineWidth',1.5,'Color',cScale(i,:), ...
        'MarkerFaceColor',cScale(i,:))
    box on;
    yyaxis right;
    plot(kappa_c,A_max(:,i).*qhat_tuned(:,i)./qref ...
        ,'--+','LineWidth',1.5,'Color',cScale(i,:))
end
set(gca,'XScale','log')
xlabel('$\kappa_\mathrm{c}$')
yyaxis left
ylabel('$A_\mathrm{max}$')
axis tight;
yyaxis right
ylabel('$(\hat{q}^\ast / \hat{q}_\mathrm{ref})_\mathrm{max}$')
savefig([savepath 'A_max.fig'])






function parsave(savepath_backup,i,j,k,qhat_tuned, qhat_tuned_std, N_sipp_tuned, qhat_mt, qhat_std_mt, N_sipp_mt)
save([savepath_backup 'qhat_tuned_' num2str(i) '_' num2str(j) '_' num2str(k) '.mat'],'qhat_tuned')
save([savepath_backup 'qhat_tuned_std_' num2str(i) '_' num2str(j) '_' num2str(k) '.mat'],'qhat_tuned_std')
save([savepath_backup 'N_sipp_tuned_' num2str(i) '_' num2str(j) '_' num2str(k) '.mat'],'N_sipp_tuned')
save([savepath_backup 'qhat_mt_' num2str(i) '_' num2str(j) '_' num2str(k) '.mat'],'qhat_mt')
save([savepath_backup 'qhat_std_mt_' num2str(i) '_' num2str(j) '_' num2str(k) '.mat'],'qhat_std_mt')
save([savepath_backup 'N_sipp_mt_' num2str(i) '_' num2str(j) '_' num2str(k) '.mat'],'N_sipp_mt')
end