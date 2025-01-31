%% Initialize variables
% Amplitude magnification factors

if sys.sigma_omega ~= 0
    % Linear case if linear mistuning enabled
    A_ref = zeros(...
        simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c, ...
        simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale, ...
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
    
    % Inverse Participation Ratio
    IPR_ref = zeros(...
        simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c, ...
        simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale, ...
              simsetup.VariationCouplingAndClearanceMCS.N_MCS);

    % Localization Factor
    LF_ref = zeros(...
        simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c, ...
        simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale, ...
              simsetup.VariationCouplingAndClearanceMCS.N_MCS);

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
% Inverse Participation Ratio
IPR = zeros(...
    simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c, ...
    simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale, ...
          simsetup.VariationCouplingAndClearanceMCS.N_MCS);

% Localization Factor
LF = zeros(...
    simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c, ...
    simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale, ...
          simsetup.VariationCouplingAndClearanceMCS.N_MCS);

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
            Resp_type_tuned(i,j) = 0; % GSR
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
                               1.01*sys_mt.r_k_mt(end),12000);
            qall =  abs(ComputeLinearResponse(r_lin,sys_mt,exc,...
                                        'mistuned','fixed_absorbers'));

            % Extract overall maximum and corresponding frequency
            [q_max,i_max] = max(max(qall,[],1),[],2);
            r_range = r_lin(i_max);

            % Localization measures
            [IPR_ref(i,j,k),LF_ref(i,j,k)] = LocalizationMeasures( ...
                qall(:,i_max),sys_mt);
            
            % Store linear amplitude magnification only for first
            % nominal clearance (only depends on coupling in linear case)
            if sys_mt.sigma_omega~=0
                A_ref(i,j,k) = q_max/sys_mt.qref;

                if j == 1
                    delta_omega_Aref(:,i,k) = sys_mt.delta_omega;
                end
            end
            
            % Removed absorbers
            % Frequency response function of max amplitude
            r_lin = linspace(0.99*sys_mt.r_k_noabs_mt(1),...
                               1.01*sys_mt.r_k_noabs_mt(end),12000);
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
            [qhat_mt_temp, qhat_mt_std_temp, N_sipp_mt_temp, IPR_temp, LF_temp] = ....
            FindResonance(sys_mt,sol,exc,r_steps,'mistuned');
            qhat_mt(i,j,k) = qhat_mt_temp;
            qhat_mt_std(i,j,k) = qhat_mt_std_temp;
            N_sipp_mt(:,i,j,k) = N_sipp_mt_temp;
            IPR(i,j,k) = IPR_temp;
            LF(i,j,k) = LF_temp;

            % If absorber malfunction: Check if impact occured in
            % malfunctioning sector
            if isfield(sys_mt,'absorber_malfunction') && ...
                 N_sipp_mt_temp(1)~=0
                
                warning('Impact detected in sector without absorber.')

            end
            
            % Save backup
            parsave(savepath_backup,i,j,k,qhat_tuned, qhat_tuned_std,...
                N_sipp_tuned, qhat_mt_temp, qhat_mt_std_temp, N_sipp_mt_temp)

            % Compute amplitude magnification factor
            A(i,j,k) = qhat_mt_temp/qhat_tuned(i,j);

            % Response type
            if isfield(sys_mt,'absorber_malfunction')

                Resp_type_mt(i,j,k) = ResponseType(N_sipp_mt_temp(2:end));
                
            else

                Resp_type_mt(i,j,k) = ResponseType(N_sipp_mt_temp);

            end
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
        [A_ref_max(i), k_max] = max(A_ref(i,1,:),[],3);
        delta_omega_Aref_max(:,i) = delta_omega_Aref(:,i,k_max);
        save([savepath 'A_ref_max.mat'],'A_ref_max')
    
        % 95% Interval
        A_ref_95(i) = quantile(squeeze(A_ref(i,1,:)),0.95);
        save([savepath 'A_ref_95.mat'],'A_ref_95')
    end

end

% Save data
save([savepath 'A.mat'],'A')
save([savepath 'A_max.mat'],'A_max')
save([savepath 'LF.mat'],'LF')
save([savepath 'IPR.mat'],'IPR')
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
    save([savepath 'LF_ref.mat'],'LF_ref')
    save([savepath 'IPR_ref.mat'],'IPR_ref')

    figure(1)
    tiledlayout(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c,...
            simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale)
    figure(2)
    tiledlayout(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c,...
            simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale)

    figure(3)
    tiledlayout(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c,...
            simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale)

    figure(4)
    tiledlayout(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c,...
            simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale)

    
    for i = 1:(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c*...
             simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale)
        
        figure(1)
        nexttile
        s=scatterhistogram(squeeze(A(i_ind(i),j_ind(i),:)), ...
            squeeze(A_ref(i_ind(i),j_ind(i),:)),...
        'HistogramDisplayStyle','bar','MarkerSize',1,'LineStyle','none');
        s.Color = {color.ies};
        s.LineWidth = 0.1;
        xlabel('$A$')
        ylabel('$A_\mathrm{ref}$')
        title(['PCC = ' num2str(round(1000*corr(squeeze(...
        A_ref(i_ind(i),j_ind(i),:)), ...
        squeeze(A(i_ind(i),j_ind(i),:))))/1000)])

        figure(2)
        nexttile
        s=scatterhistogram(squeeze(qhat_mt(i_ind(i),j_ind(i),:))/qref(i_ind(i)) ...
            ,squeeze(A_ref(i_ind(i),j_ind(i),:)),...
        'HistogramDisplayStyle','bar','MarkerSize',1,'LineStyle','none');
        s.Color = {color.ies};
        xlabel('$\mathrm{max}_j\left\{\hat{q}_j^\ast\right\} / \hat{q}_\mathrm{ref}$')
        ylabel('$A_\mathrm{ref}$')
        title(['PCC = ' num2str(round(1000*corr(squeeze(...
        A_ref(i_ind(i),j_ind(i),:)), ...
        squeeze(qhat_mt(i_ind(i),j_ind(i),:))))/1000)])

        figure(3)
        nexttile
        s=scatterhistogram(squeeze(IPR(i_ind(i),j_ind(i),:)), ...
                        squeeze(IPR_ref(i_ind(i),j_ind(i),:)),...
        'HistogramDisplayStyle','bar','MarkerSize',1,'LineStyle','none');
        s.Color = {color.ies};
        xlabel('$\mathrm{IPR}$')
        ylabel('$\mathrm{IPR}_\mathrm{ref}$')
        title(['PCC = ' num2str(round(1000*corr(squeeze(IPR(i_ind(i),j_ind(i),:)), ...
        squeeze(IPR_ref(i_ind(i),j_ind(i),:))))/1000)])


        figure(4)
        nexttile
        s=scatterhistogram(squeeze(LF(i_ind(i),j_ind(i),:)), ...
                        squeeze(LF_ref(i_ind(i),j_ind(i),:)),...
        'HistogramDisplayStyle','bar','MarkerSize',1,'LineStyle','none');
        s.Color = {color.ies};
        xlabel('$\mathrm{LF}$')
        ylabel('$\mathrm{LF}_\mathrm{ref}$')
        title(['PCC = ' num2str(round(1000*corr(squeeze(LF(i_ind(i),j_ind(i),:)), ...
        squeeze(LF_ref(i_ind(i),j_ind(i),:))))/1000)])
    end

    figure(1)
    savefig([savepath 'A_A_ref_relation.fig'])

    figure(2)
    savefig([savepath 'qhat_A_ref_relation.fig'])

    figure(3)
    savefig([savepath 'IPR_relation.fig'])

    figure(4)
    savefig([savepath 'LF_relation.fig'])

end


figure(5)
tiledlayout(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c,...
        simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale)
figure(6)
tiledlayout(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c,...
        simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale)
figure(7)
tiledlayout(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c,...
        simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale)
figure(8)
tiledlayout(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c,...
        simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale)

for i = 1:(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c*...
         simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale)
    
    [E,x] = ecdf(squeeze(A(i_ind(i),j_ind(i),:)));

    figure(5)
    nexttile;
    hold on;
    colororder({color.ies;color.analytics})
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

    figure(6)
    nexttile
    hold on;
    colororder({color.ies;color.analytics})
    yyaxis right;
    stairs(x,E,'--k','LineWidth',1.5)
    ylabel('CDF')
    yyaxis left;
    histogram(squeeze(qhat_mt(i_ind(i),j_ind(i),:))/qref(i_ind(i)),...
    'Normalization','pdf','FaceColor',color.ies,'LineStyle','none');
    box on;
    xlabel('$\mathrm{max}_j\left\{\hat{q}_j^\ast\right\} / \hat{q}_\mathrm{ref}$')
    ylabel('PDF')
    axis tight

    [E,x] = ecdf(squeeze(IPR(i_ind(i),j_ind(i),:)));

    figure(7)
    nexttile
    hold on;
    colororder({color.ies;color.analytics})
    yyaxis right;
    stairs(x,E,'--k','LineWidth',1.5)
    ylabel('CDF')
    yyaxis left;
    histogram(squeeze(IPR(i_ind(i),j_ind(i),:)),...
    'Normalization','pdf','FaceColor',color.ies,'LineStyle','none');
    box on;
    xlabel('$\mathrm{IPR}$')
    ylabel('PDF')
    axis tight

    [E,x] = ecdf(squeeze(LF(i_ind(i),j_ind(i),:)));

    figure(8)
    nexttile
    hold on;
    colororder({color.ies;color.analytics})
    yyaxis right;
    stairs(x,E,'--k','LineWidth',1.5)
    ylabel('CDF')
    yyaxis left;
    histogram(squeeze(LF(i_ind(i),j_ind(i),:)),...
    'Normalization','pdf','FaceColor',color.ies,'LineStyle','none');
    box on;
    xlabel('$\mathrm{LF}$')
    ylabel('PDF')
    axis tight

end


figure(5)
savefig([savepath 'A_PDF.fig'])
figure(6)
savefig([savepath 'qhat_PDF.fig'])
figure(7)
savefig([savepath 'IPR_PDF.fig'])
figure(8)
savefig([savepath 'LF_PDF.fig'])


% Plot response types
figure(9)
tiledlayout(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c,...
        simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale)
for i = 1:(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c*...
             simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale)

    nexttile;
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
    lstyleA = {'-o';':x';'--diamond'};
    lstyleq = {'-o';':x';'--diamond'};
    msize = [6,8,6];
else
    lstyleA = repmat({'-o';'--o';':o'}, ...
        [simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale 1]);
    lstyleq = repmat({'-o';'--o';':o'}, ...
        [simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale 1]);
    msize = 6*ones(simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale,1);
end

figure(10)
hold on;
colororder({color.analytics;color.ies})
for i = 1:simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale
    yyaxis left;
    plot(kappa_c,A_95(:,i),lstyleA{i},'LineWidth',1.5,'Color',color.analytics, ...
        'MarkerFaceColor', color.analytics, 'MarkerSize',msize(i))
    box on;
    yyaxis right;
    plot(kappa_c,(A_95(:,i).*qhat_tuned(:,i))./qref ...
        ,lstyleq{i},'LineWidth',1.5,'Color',color.ies, ...
        'MarkerFaceColor',color.ies,'HandleVisibility','off', ...
        'MarkerSize',msize(i))
end
set(gca,'XScale','log')
xlabel('$\kappa_\mathrm{c}$')
if simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale==3
    legend('$\Gamma_\mathrm{GSR}$','$\Gamma_\mathrm{OPT}$', ...
        '$\Gamma_\mathrm{SMR}$')
end
yyaxis left
ylabel('$A_{95}$')
axis tight;
yyaxis right
ylabel(['$\left[ \mathrm{max}_j \left\{ \hat{q}_j^\ast \right\} ' ...
    '/ \hat{q}_\mathrm{ref} \right]_{95}$'])
axis tight;
savefig([savepath 'A_95.fig'])

figure(11)
hold on;
colororder({color.analytics;color.ies})
for i = 1:simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale
    yyaxis left;
    plot(kappa_c,A_max(:,i),lstyleA{i},'LineWidth',1.5,'Color',color.analytics, ...
        'MarkerFaceColor',color.analytics,'MarkerSize',msize(i))
    box on;
    yyaxis right;
    plot(kappa_c,A_max(:,i).*qhat_tuned(:,i)./qref ...
        ,lstyleq{i},'LineWidth',1.5,'Color',color.ies, ...
        'MarkerFaceColor',color.ies,'HandleVisibility','off', ...
        'MarkerSize',msize(i))
end
set(gca,'XScale','log')
xlabel('$\kappa_\mathrm{c}$')
if simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale==3
    legend('$\Gamma_\mathrm{GSR}$','$\Gamma_\mathrm{OPT}$', ...
        '$\Gamma_\mathrm{SMR}$')
end
yyaxis left
ylabel('$A_\mathrm{max}$')
axis tight;
yyaxis right
ylabel(['$\mathrm{max}_\mathrm{MCS} \left\{ \mathrm{max}_j \left\{ \hat{q}_j^\ast \right\} ' ...
    '/ \hat{q}_\mathrm{ref} \right\}$'])
axis tight;
savefig([savepath 'A_max.fig'])

figure(13)
tiledlayout(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c,...
        simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale)

figure(14)
tiledlayout(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c,...
        simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale)

figure(15)
tiledlayout(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c,...
        simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale)

figure(16)
tiledlayout(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c,...
        simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale)

for i = 1:(simsetup.VariationCouplingAndClearanceMCS.Number_kappa_c*...
             simsetup.VariationCouplingAndClearanceMCS.Number_GammaScale)
        
    figure(13)
    nexttile
    s=scatterhistogram(squeeze(A(i_ind(i),j_ind(i),:)), ...
        squeeze(IPR(i_ind(i),j_ind(i),:)),...
    'HistogramDisplayStyle','bar','MarkerSize',1,'LineStyle','none');
    s.Color = {color.ies};
    s.LineWidth = 0.1;
    xlabel('$A$')
    ylabel('$\mathrm{IPR}$')
    title(['PCC = ' num2str(round(1000*corr(squeeze( ...
        IPR(i_ind(i),j_ind(i),:)), ...
    squeeze(A(i_ind(i),j_ind(i),:))))/1000)])

    figure(14)
    nexttile
    s=scatterhistogram(squeeze(A(i_ind(i),j_ind(i),:)), ...
        squeeze(LF(i_ind(i),j_ind(i),:)),...
    'HistogramDisplayStyle','bar','MarkerSize',1,'LineStyle','none');
    s.Color = {color.ies};
    s.LineWidth = 0.1;
    xlabel('$A$')
    ylabel('$\mathrm{LF}$')
    title(['PCC = ' num2str(round(1000*corr(squeeze( ...
        LF(i_ind(i),j_ind(i),:)), ...
    squeeze(A(i_ind(i),j_ind(i),:))))/1000)])

    figure(15)
    nexttile
    s=scatterhistogram(squeeze(qhat_mt(i_ind(i),j_ind(i),:))/qref(i_ind(i)), ...
        squeeze(IPR(i_ind(i),j_ind(i),:)),...
    'HistogramDisplayStyle','bar','MarkerSize',1,'LineStyle','none');
    s.Color = {color.ies};
    s.LineWidth = 0.1;
    xlabel('$\mathrm{max}_j\left\{\hat{q}_j^\ast\right\} / \hat{q}_\mathrm{ref}$')
    ylabel('$\mathrm{IPR}$')
    title(['PCC = ' num2str(round(1000*corr(squeeze( ...
        IPR(i_ind(i),j_ind(i),:)), ...
    squeeze(qhat_mt(i_ind(i),j_ind(i),:))/qref(i_ind(i))))/1000)])

    figure(16)
    nexttile
    s=scatterhistogram(squeeze(qhat_mt(i_ind(i),j_ind(i),:))/qref(i_ind(i)), ...
        squeeze(LF(i_ind(i),j_ind(i),:)),...
    'HistogramDisplayStyle','bar','MarkerSize',1,'LineStyle','none');
    s.Color = {color.ies};
    s.LineWidth = 0.1;
    xlabel('$\mathrm{max}_j\left\{\hat{q}_j^\ast\right\} / \hat{q}_\mathrm{ref}$')
    ylabel('$\mathrm{LF}$')
    title(['PCC = ' num2str(round(1000*corr(squeeze( ...
        LF(i_ind(i),j_ind(i),:)), ...
    squeeze(qhat_mt(i_ind(i),j_ind(i),:))/qref(i_ind(i))))/1000)])

        
end
figure(13)
savefig([savepath 'A_IPR_relation.fig'])
figure(14)
savefig([savepath 'A_LF_relation.fig'])
figure(15)
savefig([savepath 'qhat_IPR_relation.fig'])
figure(16)
savefig([savepath 'qhat_LF_relation.fig'])



function parsave(savepath_backup,i,j,k,qhat_tuned, qhat_tuned_std, N_sipp_tuned, qhat_mt, qhat_std_mt, N_sipp_mt)
save([savepath_backup 'qhat_tuned_' num2str(i) '_' num2str(j) '_' num2str(k) '.mat'],'qhat_tuned')
save([savepath_backup 'qhat_tuned_std_' num2str(i) '_' num2str(j) '_' num2str(k) '.mat'],'qhat_tuned_std')
save([savepath_backup 'N_sipp_tuned_' num2str(i) '_' num2str(j) '_' num2str(k) '.mat'],'N_sipp_tuned')
save([savepath_backup 'qhat_mt_' num2str(i) '_' num2str(j) '_' num2str(k) '.mat'],'qhat_mt')
save([savepath_backup 'qhat_std_mt_' num2str(i) '_' num2str(j) '_' num2str(k) '.mat'],'qhat_std_mt')
save([savepath_backup 'N_sipp_mt_' num2str(i) '_' num2str(j) '_' num2str(k) '.mat'],'N_sipp_mt')
end