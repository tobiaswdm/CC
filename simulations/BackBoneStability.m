%% Initialize variables
qhat_practically_stable = zeros(simsetup.BackBoneStability.Nkappac, ...
    simsetup.BackBoneStability.Nxi);
qhat_stable = zeros(simsetup.BackBoneStability.Nkappac, ...
    simsetup.BackBoneStability.Nxi);
qhat_unstable = zeros(simsetup.BackBoneStability.Nkappac, ...
    simsetup.BackBoneStability.Nxi);
qhat_unstable_synchloss = zeros(simsetup.BackBoneStability.Nkappac, ...
    simsetup.BackBoneStability.Nxi);
qhat_unstable_modulation = zeros(simsetup.BackBoneStability.Nkappac, ...
    simsetup.BackBoneStability.Nxi);
lambda_real = zeros(simsetup.BackBoneStability.Nkappac, ...
    simsetup.BackBoneStability.Nxi);
lambda_imag = zeros(simsetup.BackBoneStability.Nkappac, ...
    simsetup.BackBoneStability.Nxi);

% Samples of coupling strength
kappa_c = logspace(log10(simsetup.BackBoneStability.kappac_range(1)), ...
                  log10(simsetup.BackBoneStability.kappac_range(2)), ...
                  simsetup.BackBoneStability.Nkappac);


%% Determine Backbone

% Initail estimation of backbone
[Gamma_scale,q_scale,Gamma_opt] = TunedBackbone(sys,'stable');

% Sample clearance-nomralized amplitude
xi = logspace(log10(min(q_scale)/Gamma_opt),log10(5.1),...
              simsetup.BackBoneStability.Nxi);

% Auxiliary variable
rho = (2/pi)*(1-sys.eN)/(sys.eN+1);
% Backbone
theta = (1+sqrt((1+rho^2)*xi.^2 - rho^2))/(1+rho^2);
varpi_bb = sqrt(xi./((1-sys.epsilon_a)*xi + ...
    8*sys.epsilon_a*theta.*((theta-1)./xi)/pi^2));
Delta = acos((theta-1)./xi);
Gamma_scale_bb = 2*sys.D./abs((-(1-sys.epsilon_a)*varpi_bb.^2 + ...
    2*sys.D*1i*varpi_bb + 1).*xi - ...
    8*sys.epsilon_a*theta.*varpi_bb.^2 .* exp(-1i*Delta)/pi^2);

%% Study Stability along Backbone

for i = 1:simsetup.BackBoneStability.Nkappac

    disp(['Inter-sector coupling strength ' num2str(i) ...
            ' of ' num2str(simsetup.BackBoneStability.Nkappac)])

    % Assign coupling
    sys.kappa_c = kappa_c(i);
    
    parfor j = 1:simsetup.BackBoneStability.Nxi

        sys_temp = sys;

        sys_temp.Gamma_Scale = Gamma_scale_bb(j);

        [sys_temp,exc_temp] = BuildSystem(sys_temp,exc,'tuned');
        
        % Put in pseudo contour array
        c = [1 , varpi_bb(j)*sys_temp.r_k(exc.k+1);
             1,  xi(j)];
        
        % Stability Analysis
        [qhat_practically_stable(i,j),qhat_stable(i,j),qhat_unstable(i,j),~, ...
        qhat_unstable_synchloss(i,j), qhat_unstable_modulation(i,j), ...
        lambda] = ...
        StabilityAnalysisGSR(c,sys_temp,sol,exc_temp);
        
        % Real and imaginary part of leading eigenvalue
        lambda_real(i,j) = real(lambda);
        lambda_imag(i,j) = imag(lambda);
    end

    kc_plot = sys.kappa_c*ones(1,simsetup.BackBoneStability.Nxi);
    
    if i == 1
        figure(1)
        hold on; box on;
        scatter(kc_plot(~isnan(qhat_unstable(i,:))), ...
            Gamma_scale_bb(~isnan(qhat_unstable(i,:))),20, ...
            'MarkerFaceColor',color.show,...
            'MarkerEdgeColor','k','Displayname','Unstable')
        scatter(kc_plot(~isnan(qhat_stable(i,:))), ...
            Gamma_scale_bb(~isnan(qhat_stable(i,:))), ...
            20,'MarkerFaceColor',myColors('cyan'),...
            'MarkerEdgeColor','k','Displayname','L. A. Stable')
        scatter(kc_plot(~isnan(qhat_practically_stable(i,:))), ...
            Gamma_scale_bb(~isnan(qhat_practically_stable(i,:))), ...
            40,'pentagram','MarkerFaceColor',myColors('green'), ...
            'MarkerEdgeColor','k','Displayname','Pract. Stable')
        set(gca,'YScale','log','XScale','log')
        xlabel('$\kappa_\mathrm{c}$')
        ylabel('$\Gamma / \hat{q}_\mathrm{ref}$')

        figure(2)
        hold on; box on;
        scatter(kc_plot(~isnan(qhat_unstable(i,:))), ...
            Gamma_scale_bb(~isnan(qhat_unstable(i,:))),20, ...
            'MarkerFaceColor',color.show,...
            'MarkerEdgeColor','k','Displayname','Unstable', ...
            'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25)
        scatter(kc_plot(~isnan(qhat_stable(i,:)) &...
                isnan(qhat_unstable_synchloss(i,:))), ...
                Gamma_scale_bb(~isnan(qhat_stable(i,:)) &...
                isnan(qhat_unstable_synchloss(i,:))), ...
                20,'MarkerFaceColor',myColors('cyan'),...
                'MarkerEdgeColor','k','Displayname','L. A. Stable', ...
                'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25)
        scatter(kc_plot(~isnan(qhat_unstable_synchloss(i,:))), ...
                Gamma_scale_bb(~isnan(qhat_unstable_synchloss(i,:))), ...
                20,'MarkerFaceColor',myColors('magenta'),...
                'MarkerEdgeColor','k','Displayname','Criterion 1')
        scatter(kc_plot(~isnan(qhat_practically_stable(i,:))), ...
            Gamma_scale_bb(~isnan(qhat_practically_stable(i,:))), ...
            40,'pentagram','MarkerFaceColor',myColors('green'), ...
            'MarkerEdgeColor','k','Displayname','Pract. Stable', ...
            'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25)
        set(gca,'YScale','log','XScale','log')
        xlabel('$\kappa_\mathrm{c}$')
        ylabel('$\Gamma / \hat{q}_\mathrm{ref}$')
        title('Criterion 1')

        figure(3)
        hold on; box on;
        scatter(kc_plot(~isnan(qhat_unstable(i,:))), ...
            Gamma_scale_bb(~isnan(qhat_unstable(i,:))),20, ...
            'MarkerFaceColor',color.show,...
            'MarkerEdgeColor','k','Displayname','Unstable', ...
            'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25)
        scatter(kc_plot(~isnan(qhat_stable(i,:)) &...
                isnan(qhat_unstable_modulation(i,:))), ...
                Gamma_scale_bb(~isnan(qhat_stable(i,:))&...
                isnan(qhat_unstable_modulation(i,:))), ...
                20,'MarkerFaceColor',myColors('cyan'),...
                'MarkerEdgeColor','k','Displayname','L. A. Stable', ...
                'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25)
        scatter(kc_plot(~isnan(qhat_unstable_modulation(i,:))), ...
                Gamma_scale_bb(~isnan(qhat_unstable_modulation(i,:))), ...
                20,'MarkerFaceColor',myColors('magenta'),...
                'MarkerEdgeColor','k','Displayname','Criterion 2')
        scatter(kc_plot(~isnan(qhat_practically_stable(i,:))), ...
            Gamma_scale_bb(~isnan(qhat_practically_stable(i,:))), ...
            40,'pentagram','MarkerFaceColor',myColors('green'), ...
            'MarkerEdgeColor','k','Displayname','Pract. Stable', ...
            'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25)
        set(gca,'YScale','log','XScale','log')
        xlabel('$\kappa_\mathrm{c}$')
        ylabel('$\Gamma / \hat{q}_\mathrm{ref}$')
        title('Criterion 2')

        figure(4)
        hold on; box on;
        scatter(kc_plot(~isnan(qhat_unstable(i,:))), ...
            Gamma_scale_bb(~isnan(qhat_unstable(i,:))),20, ...
            'MarkerFaceColor',color.show,...
            'MarkerEdgeColor','k','Displayname','Unstable', ...
            'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25)
        scatter(kc_plot(~isnan(qhat_stable(i,:)) & ...
                ~isnan(qhat_unstable_synchloss(i,:))), ...
                Gamma_scale_bb(~isnan(qhat_stable(i,:)) &...
                ~isnan(qhat_unstable_synchloss(i,:))), ...
                20,'MarkerFaceColor',myColors('cyan'),...
                'MarkerEdgeColor','k','Displayname','L. A. Stable', ...
                'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25)
        scatter(kc_plot(~isnan(qhat_unstable_modulation(i,:)) & ...
                         isnan(qhat_unstable_synchloss(i,:))), ...
                Gamma_scale_bb(~isnan(qhat_unstable_modulation(i,:)) & ...
                         isnan(qhat_unstable_synchloss(i,:))), ...
                20,'MarkerFaceColor',myColors('magenta'),...
                'MarkerEdgeColor','k','Displayname','Modulation')
        scatter(kc_plot(~isnan(qhat_practically_stable(i,:))), ...
            Gamma_scale_bb(~isnan(qhat_practically_stable(i,:))), ...
            40,'pentagram','MarkerFaceColor',myColors('green'), ...
            'MarkerEdgeColor','k','Displayname','Pract. Stable', ...
            'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25)
        set(gca,'YScale','log','XScale','log')
        xlabel('$\kappa_\mathrm{c}$')
        ylabel('$\Gamma / \hat{q}_\mathrm{ref}$')
        title('Breather Candidates')

        figure(5)
        hold on; box on;
        scatter(kc_plot(~isnan(qhat_unstable(i,:))), ...
            Gamma_scale_bb(~isnan(qhat_unstable(i,:))),20, ...
            'MarkerFaceColor',color.show,...
            'MarkerEdgeColor','k','Displayname','Unstable')
        scatter(kc_plot(~isnan(qhat_stable(i,:))), ...
            Gamma_scale_bb(~isnan(qhat_stable(i,:))), ...
            20,'MarkerFaceColor',myColors('cyan'),...
            'MarkerEdgeColor','k','Displayname','L. A. Stable')
        scatter(kc_plot(~isnan(qhat_practically_stable(i,:))), ...
            Gamma_scale_bb(~isnan(qhat_practically_stable(i,:))), ...
            20,'MarkerFaceColor',myColors('cyan'), ...
            'MarkerEdgeColor','k','HandleVisibility','off')
        set(gca,'YScale','log','XScale','log')
        xlabel('$\kappa_\mathrm{c}$')
        ylabel('$\Gamma / \hat{q}_\mathrm{ref}$')

       

    else
        figure(1)
        hold on;
        scatter(kc_plot(~isnan(qhat_unstable(i,:))), ...
            Gamma_scale_bb(~isnan(qhat_unstable(i,:))),20, ...
            'MarkerFaceColor',color.show,...
            'MarkerEdgeColor','k','HandleVisibility','off')
        scatter(kc_plot(~isnan(qhat_stable(i,:))), ...
            Gamma_scale_bb(~isnan(qhat_stable(i,:))), ...
            20,'MarkerFaceColor',myColors('cyan'),...
            'MarkerEdgeColor','k','HandleVisibility','off')
        scatter(kc_plot(~isnan(qhat_practically_stable(i,:))), ...
            Gamma_scale_bb(~isnan(qhat_practically_stable(i,:))), ...
            40,'pentagram','MarkerFaceColor',myColors('green'), ...
            'MarkerEdgeColor','k','HandleVisibility','off')

        figure(2)
        hold on; box on;
        scatter(kc_plot(~isnan(qhat_unstable(i,:))), ...
            Gamma_scale_bb(~isnan(qhat_unstable(i,:))),20, ...
            'MarkerFaceColor',color.show,...
            'MarkerEdgeColor','k','HandleVisibility','off', ...
            'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25)
        scatter(kc_plot(~isnan(qhat_stable(i,:)) &...
                isnan(qhat_unstable_synchloss(i,:))), ...
                Gamma_scale_bb(~isnan(qhat_stable(i,:)) &...
                isnan(qhat_unstable_synchloss(i,:))), ...
                20,'MarkerFaceColor',myColors('cyan'),...
                'MarkerEdgeColor','k','HandleVisibility','off', ...
                'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25)
        scatter(kc_plot(~isnan(qhat_unstable_synchloss(i,:))), ...
                Gamma_scale_bb(~isnan(qhat_unstable_synchloss(i,:))), ...
                20,'MarkerFaceColor',myColors('magenta'),...
                'MarkerEdgeColor','k','HandleVisibility','off')
        scatter(kc_plot(~isnan(qhat_practically_stable(i,:))), ...
            Gamma_scale_bb(~isnan(qhat_practically_stable(i,:))), ...
            40,'pentagram','MarkerFaceColor',myColors('green'), ...
            'MarkerEdgeColor','k','HandleVisibility','off', ...
            'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25)

        figure(3)
        hold on; box on;
        scatter(kc_plot(~isnan(qhat_unstable(i,:))), ...
            Gamma_scale_bb(~isnan(qhat_unstable(i,:))),20, ...
            'MarkerFaceColor',color.show,...
            'MarkerEdgeColor','k','HandleVisibility','off', ...
            'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25)
        scatter(kc_plot(~isnan(qhat_stable(i,:)) &...
                isnan(qhat_unstable_modulation(i,:))), ...
                Gamma_scale_bb(~isnan(qhat_stable(i,:))&...
                isnan(qhat_unstable_modulation(i,:))), ...
                20,'MarkerFaceColor',myColors('cyan'),...
                'MarkerEdgeColor','k','HandleVisibility','off', ...
                'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25)
        scatter(kc_plot(~isnan(qhat_unstable_modulation(i,:))), ...
                Gamma_scale_bb(~isnan(qhat_unstable_modulation(i,:))), ...
                20,'MarkerFaceColor',myColors('magenta'),...
                'MarkerEdgeColor','k','HandleVisibility','off')
        scatter(kc_plot(~isnan(qhat_practically_stable(i,:))), ...
            Gamma_scale_bb(~isnan(qhat_practically_stable(i,:))), ...
            40,'pentagram','MarkerFaceColor',myColors('green'), ...
            'MarkerEdgeColor','k','HandleVisibility','off', ...
            'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25)

        figure(4)
        hold on; box on;
        scatter(kc_plot(~isnan(qhat_unstable(i,:))), ...
            Gamma_scale_bb(~isnan(qhat_unstable(i,:))),20, ...
            'MarkerFaceColor',color.show,...
            'MarkerEdgeColor','k','HandleVisibility','off', ...
            'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25)
        scatter(kc_plot(~isnan(qhat_stable(i,:)) &...
                ~isnan(qhat_unstable_synchloss(i,:))), ...
                Gamma_scale_bb(~isnan(qhat_stable(i,:)) &...
                ~isnan(qhat_unstable_synchloss(i,:))), ...
                20,'MarkerFaceColor',myColors('cyan'),...
                'MarkerEdgeColor','k','HandleVisibility','off', ...
                'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25)
        scatter(kc_plot(~isnan(qhat_unstable_modulation(i,:)) & ...
                         isnan(qhat_unstable_synchloss(i,:))), ...
                Gamma_scale_bb(~isnan(qhat_unstable_modulation(i,:)) & ...
                         isnan(qhat_unstable_synchloss(i,:))), ...
                20,'MarkerFaceColor',myColors('magenta'),...
                'MarkerEdgeColor','k','HandleVisibility','off')
        scatter(kc_plot(~isnan(qhat_practically_stable(i,:))), ...
            Gamma_scale_bb(~isnan(qhat_practically_stable(i,:))), ...
            40,'pentagram','MarkerFaceColor',myColors('green'), ...
            'MarkerEdgeColor','k','HandleVisibility','off', ...
            'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25)
        set(gca,'YScale','log','XScale','log')
        xlabel('$\kappa_\mathrm{c}$')
        ylabel('$\Gamma / \hat{q}_\mathrm{ref}$')
        title('Breather Candidates')

        figure(5)
        hold on;
        scatter(kc_plot(~isnan(qhat_unstable(i,:))), ...
            Gamma_scale_bb(~isnan(qhat_unstable(i,:))),20, ...
            'MarkerFaceColor',color.show,...
            'MarkerEdgeColor','k','HandleVisibility','off')
        scatter(kc_plot(~isnan(qhat_stable(i,:))), ...
            Gamma_scale_bb(~isnan(qhat_stable(i,:))), ...
            20,'MarkerFaceColor',myColors('cyan'),...
            'MarkerEdgeColor','k','HandleVisibility','off')
        scatter(kc_plot(~isnan(qhat_practically_stable(i,:))), ...
            Gamma_scale_bb(~isnan(qhat_practically_stable(i,:))), ...
            20,'MarkerFaceColor',myColors('cyan'), ...
            'MarkerEdgeColor','k','HandleVisibility','off')

    end


end

[KK,GG] = meshgrid(kappa_c,Gamma_scale_bb);

% Leading Ev
figure(6);
tiledlayout(1,2)
ax1=nexttile;
hold on;
surf(KK,GG,lambda_real','LineStyle','none')
xlabel('$\kappa_\mathrm{c}$')
ylabel('$\Gamma / \hat{q}_\mathrm{ref}$')
title('Real of leading EV')
colormap(ax1,redblue)
xlabel('$\kappa_\mathrm{c}$')
ylabel('$\Gamma/\hat{q}_\mathrm{ref}$')
h1=colorbar(ax1);
h1.Label.Interpreter = 'latex';
h1.Label.String = "Re - Leading EV";
clim(ax1,max(abs(lambda_real),[],'all')*[-1 1])
set(gca,'YDir','normal')
set(gca,'XScale','log')
set(gca,'YScale','log')
axis tight;
view([0 90])
hold off
ax2=nexttile;
hold on;
surf(KK,GG,lambda_imag','LineStyle','none')
xlabel('$\kappa_\mathrm{c}$')
ylabel('$\Gamma / \hat{q}_\mathrm{ref}$')
title('Imag of leading EV')
colormap(ax2,plasma)
xlabel('$\kappa_\mathrm{c}$')
ylabel('$\Gamma/\hat{q}_\mathrm{ref}$')
h2=colorbar(ax2);
h2.Label.Interpreter = 'latex';
h2.Label.String = "Im - Leading EV";
%clim(ax2,[min(abs(lambda_imag),[],'all'),max(abs(lambda_imag),[],'all')])
set(ax2,'ColorScale','linear')
set(gca,'YDir','normal')
set(gca,'XScale','log')
set(gca,'YScale','log')
axis tight;
view([0 90])
hold off;
savefig([savepath 'leading_EV.fig'])

%% Save data 

% Figures
figure(1)
savefig([savepath 'backbone_stability.fig'])
figure(2)
savefig([savepath 'synchloss.fig'])
figure(3)
savefig([savepath 'modulation.fig'])
figure(4)
savefig([savepath 'breather.fig'])
figure(5)
savefig([savepath 'backbone_pure_stability.fig'])

% Data
save([savepath 'xi.mat'],'xi')
save([savepath 'varpi_bb.mat'],'varpi_bb')
save([savepath 'Gamma_scale_bb.mat'],'Gamma_scale_bb')
save([savepath 'kappa_c.mat'],'kappa_c')
save([savepath 'qhat_stable.mat'],'qhat_stable')
save([savepath 'qhat_unstable.mat'],'qhat_unstable')
save([savepath 'qhat_practically_stable.mat'],'qhat_practically_stable')
save([savepath 'qhat_unstable_synchloss.mat'],'qhat_unstable_synchloss')
save([savepath 'qhat_unstable_modulation.mat'],'qhat_unstable_modulation')