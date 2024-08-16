% Auxilliary Variable
rho = (2/pi) * (1-sys.eN) / (1+sys.eN);
% Minimum amplitude - turning point of SIM
% Including safety of 0.1% higher amplitude
xi_min = 1.001 * rho/sqrt(1+rho^2);

% Clearance normalized amplitudes
xi = logspace(log10(xi_min),...
    log10(simsetup.FollowBifurcationFRSofLSR.xi_max),...
    simsetup.FollowBifurcationFRSofLSR.Nxi);
% Excitation frequencies
r = linspace(simsetup.FollowBifurcationFRSofLSR.r_range(1),...
    simsetup.FollowBifurcationFRSofLSR.r_range(2), ...
    simsetup.FollowBifurcationFRSofLSR.Nr);
% Bifurcation Point of FRS
GammaScaleBif = zeros(1,simsetup.FollowBifurcationFRSofLSR.Nkappac);

% Samples of coupling strength
kappa_c = logspace(log10(simsetup.FollowBifurcationFRSofLSR.kappac_range(1)), ...
                  log10(simsetup.FollowBifurcationFRSofLSR.kappac_range(2)), ...
                  simsetup.FollowBifurcationFRSofLSR.Nkappac);


for i=1:simsetup.FollowBifurcationFRSofLSR.Nkappac
    disp(['Coupling stiffness ' num2str(i) ' of ' ...
        num2str(simsetup.FollowBifurcationFRSofLSR.Nkappac)])
    % Assign coupling
    sys.kappa_c = kappa_c(i);
    [sys,exc] = BuildSystem(sys,exc,'tuned');
    
    % FRS
    [Gamma_Scale,~,~] = SingleSectorFRS(xi,r,sys,exc,'tuned');
    
    % Extract the value of the local maximum that forms the first nonlinear
    % resonance branch
    TF = imregionalmax(Gamma_Scale,8);
    TFvec = sum(TF(2:end-1,2:end-1),2)~=0;
    j = find(TFvec)+1;
    TF = TF(j(1),:);
    Gamma_Scale = Gamma_Scale(j(1),:);
    GammaScaleBif(i) = Gamma_Scale(TF);
    
end

% Get optimum clearance
[~,~,Gamma_opt] = TunedBackbone(sys,'stable');


figure(1);
yline(Gamma_opt,'--k','LineWidth',1.5)
hold on;
plot(kappa_c,GammaScaleBif,'LineWidth',1.5,'Color',color.ies)
box on;
axis tight;
xlabel('$\kappa_\mathrm{c}$')
ylabel('$\Gamma / \hat{q}_\mathrm{ref}$')
set(gca,'XScale','log')
set(gca,'YScale','log')
savefig([savepath 'GammaScaleBif.fig'])
