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
    log10(simsetup.LocalizationSingleSectorAnalytical.xi_max),...
    simsetup.LocalizationSingleSectorAnalytical.Nxi);
r = linspace(simsetup.LocalizationSingleSectorAnalytical.r_range(1),...
    simsetup.LocalizationSingleSectorAnalytical.r_range(2), ...
    simsetup.LocalizationSingleSectorAnalytical.Nr);

% ESIM
[Gamma_Scale,~,Xi,R] = SingleSectorESIM(xi,r,sys,exc,'tuned');

% Plot
figure(2);
surf(R,Xi,Gamma_Scale,'EdgeAlpha',0)
hold on;
title('Tuned System')
box on;
colormap turbo
xlabel('$r$')
ylabel('$\xi$')
zlabel('$\Gamma/\hat{q}_\mathrm{ref}$')
set(gca,'YScale','log')
axis tight;

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

%% Compute example of mistuned manifold

[sys_mt,exc_mt] = BuildSystem(sys,exc,'mistuned');
[Gamma_Scale_mt,~,Xi,R] = SingleSectorESIM(xi,r,sys_mt,exc_mt,'mistuned');

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

%% Different realizations of FRFs at fixed clearances
xi_max = 0;
for i = 1:simsetup.LocalizationSingleSectorAnalytical.N_MCS
    
    [sys_mt,exc_mt] = BuildSystem(sys,exc,'mistuned');
    [Gamma_Scale_mt,~,Xi,R] = SingleSectorESIM(xi,r,sys_mt,exc_mt,'mistuned');

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
    xi_max_temp = max(xi_max,max(c(2,c(2,:)~=ceil(c(2,:)))));
    if xi_max_temp > xi_max
        xi_max = xi_max_temp;
    end
    

end

figure(6)
hold on;
c=contour(R,Xi,Gamma_Scale,sys.Gamma_Scale*[1 1],...
        'LineWidth',1.5,'EdgeColor',color.ies,...
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