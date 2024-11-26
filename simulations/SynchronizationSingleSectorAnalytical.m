%% Build tuned system

% Tuned
[sys,exc] = BuildSystem(sys,exc,'tuned');

% Draw Dispersion Diagram
DrawDispersion(sys,color,savepath);

%% Compute tuned manifold

% Auxilliary Variable
rho = (2/pi) * (1-sys.eN) / (1+sys.eN);
% Minimum amplitude - turning point of SIM
% Including safety of 0.1% higher amplitude
xi_min = 1.001 * rho/sqrt(1+rho^2);

% Clearance normalized amplitudes
xi = logspace(log10(xi_min),...
    log10(simsetup.SynchronizationSingleSectorAnalytical.xi_max),...
    simsetup.SynchronizationSingleSectorAnalytical.Nxi);
r = linspace(simsetup.SynchronizationSingleSectorAnalytical.r_range(1),...
    simsetup.SynchronizationSingleSectorAnalytical.r_range(2), ...
    simsetup.SynchronizationSingleSectorAnalytical.Nr);

% FRS
[Gamma_Scale,Xi,R] = SingleSectorFRS(xi,r,sys,exc,'tuned');

% Plot
figure(2);
surf(R,Xi,Gamma_Scale,'LineStyle','none')
hold on; box on;
contour3(R,Xi,Gamma_Scale,4,'LineWidth',1,'Color',color.analytics)
for i = 1:length(sys.r_k)
    xline(sys.r_k(i),'--k','LineWidth',1.5,'Color',[1 1 1])
end
hold off;
title('Tuned System')
colormap(jet)
xlabel('$r$')
ylabel('$\xi$')
zlabel('$\Gamma/\hat{q}_\mathrm{ref}$')
h=colorbar;
h.Label.Interpreter = 'latex';
h.Label.String = "$\Gamma/\hat{q}_\mathrm{ref}$";
set(gca,'YScale','log')
axis tight;
view([0 90])

figure(3);
for i = 1:length(sys.r_k)
    xline(sys.r_k(i),'--k',['$r_' num2str(i-1) '$'],'LineWidth',1.5, ...
        'Interpreter','Latex','LabelOrientation','horizontal')
end
hold on;
contour(R,Xi,Gamma_Scale,10,'LineWidth',1.5)
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

%% Compute example of mistuned manifold

[sys_mt,exc_mt] = BuildSystem(sys,exc,'mistuned');
[Gamma_Scale_mt,Xi,R] = SingleSectorFRS(xi,r,sys_mt,exc_mt,'mistuned');

% Plot
figure(4);
surf(R,Xi,Gamma_Scale_mt,'EdgeAlpha',0)
hold on;
title('Example mistuned System')
box on;
colormap turbo
xlabel('$r$')
ylabel('$\xi$')
zlabel('$\Gamma/\hat{q}_\mathrm{ref}$')
set(gca,'YScale','log')
axis tight;

figure(5);
contour(R,Xi,Gamma_Scale_mt,10,'LineWidth',1.5)
hold on;
h=colorbar;
h.Label.Interpreter = 'latex';
h.Label.String = "$\Gamma/\hat{q}_\mathrm{ref}$";
title('Example mistuned System')
box on;
colormap turbo
xlabel('$r$')
ylabel('$\xi$')
set(gca,'YScale','log')
axis tight;

% Linear FRF
q_fixed = ComputeLinearResponse(r,sys,exc,'tuned','fixed_absorbers');
q_fixed = abs(q_fixed(1,:));
q_removed = ComputeLinearResponse(r,sys,exc,'tuned','removed_absorbers');
q_removed = abs(q_removed(1,:));

figure(7);
hold on;
plot(r,q_fixed/sys.qref,...
    'LineWidth',.5,'Color',color.reference,'DisplayName', ...
    'Fixed abs.')
plot(r,q_removed/sys.qref,'-.',...
    'LineWidth',.5,'Color',color.reference,'DisplayName', ...
    'Removed abs.')
axis tight;


%% Different realizations of FRFs at fixed clearances
xi_max = 0;
for i = 1:simsetup.SynchronizationSingleSectorAnalytical.N_MCS
    
    [sys_mt,exc_mt] = BuildSystem(sys,exc,'mistuned');
    [Gamma_Scale_mt,Xi,R] = SingleSectorFRS(xi,r,sys_mt,exc_mt,'mistuned');

    figure(6)
    hold on;
    if i == 1
        c = contour(R,Xi,Gamma_Scale_mt,sys.Gamma_Scale*[1 1],...
        'LineWidth',.5,'EdgeColor',color.background,....
        'DisplayName','Mistuned');
    else
        c = contour(R,Xi,Gamma_Scale_mt,sys.Gamma_Scale*[1 1],...
        'LineWidth',.5,'EdgeColor',color.background,....
        'HandleVisibility','off');
    end

    % Recover maximum amplitude
    [qhat_max,~,~,r_plot] = ...
        LocalizedFrequencyAmplitudeCurve(c,sys_mt,exc_mt,'single', ...
        'mistuned');

    figure(7);
    hold on;
    if i == 1
        plot(r_plot,qhat_max/sys.qref,...
            'LineWidth',.5,'Color',color.background,'DisplayName', ...
            'Mistuned')
    else
        plot(r_plot,qhat_max/sys.qref,...
            'LineWidth',.5,'Color',color.background,...
            'HandleVisibility', 'off')
    end
    axis tight;
    

end

figure(6)
hold on;
c=contour(R,Xi,Gamma_Scale,sys.Gamma_Scale*[1 1],...
        'LineWidth',2.5,'EdgeColor',color.ies,...
        'DisplayName','Tuned');
xi_max_temp = max(xi_max,max(c(2,c(2,:)~=ceil(c(2,:)))));
if xi_max_temp > xi_max
    xi_max = xi_max_temp;
end
axis tight;
set(gca,'YScale','log')
box on
xlabel('$r$')
ylabel('$\xi$')
ylim([xi_min,xi_max])
legend;
title(['$\Gamma / \hat{q}_\mathrm{ref} = ' num2str(sys.Gamma_Scale) '$'])

% Recover maximum amplitude
[qhat_max,qhat_max_violated,~,r_plot] = ...
    LocalizedFrequencyAmplitudeCurve(c,sys,exc,'single','tuned');

figure(7)
hold on;
plot(r_plot,qhat_max/sys.qref,...
            'LineWidth',2.5,'Color',color.ies,'DisplayName', ...
            'Tuned')
plot(r_plot,qhat_max_violated/sys.qref,':',...
            'LineWidth',2.5,'Color',color.show,'DisplayName', ...
            'Tuned - Viol. kin. constr.')
axis tight;
set(gca,'YScale','log')
box on
xlabel('$r$')
ylabel('$\hat{q}/\hat{q}_\mathrm{ref}$')
legend;
title(['$\Gamma / \hat{q}_\mathrm{ref} = ' num2str(sys.Gamma_Scale) '$'])