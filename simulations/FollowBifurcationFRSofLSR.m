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
% Maximum Amplitude
MaximumAmplitude = nan(1,simsetup.FollowBifurcationFRSofLSR.Nkappac);

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
    [FRS,Xi,~] = SingleSectorFRS(xi,r,sys,exc,'tuned');
    
    % Extract the value of the local maximum that forms the first nonlinear
    % resonance branch
    TF = imregionalmax(FRS,8);
    TFvec = sum(TF(2:end-1,2:end-1),2)~=0;
    j = find(TFvec)+1;
    TF = TF(j(1),:);
    Gamma_Scale = FRS(j(1),:);
    GammaScaleBif(i) = Gamma_Scale(TF);
    
    % Check if Clearance is lower than bifurcation point
    if GammaScaleBif(i) >= sys.Gamma_Scale
        % Determine frequency Amplitude Curve
        % Get Level curves at clearance
        c = contourc(r,xi,FRS',[sys.Gamma_Scale sys.Gamma_Scale]);
        
        % Determine max amplitude FRF
        [qhat_max,~,qhat_syn,r_ana] = ...
            LocalizedFrequencyAmplitudeCurve(c,sys,exc,'single','tuned');
              
        % Does any solution fulfil kinematic constraint?
        if any(~isnan(qhat_max))
            
            % Find all boundaries between branches
            if any(isnan(qhat_max))
                
                % All indices
                indices = 1:length(qhat_max);

                j = true(size(indices));
                
                % Determine indices of all kinematically accessible 
                % solutions
                ind = find(~isnan(qhat_max));
                % Determine the kineamtically accessible solution with
                % minimum frequency
                [~,i_min] = min(r_ana(ind));
                % Index in full solution branch(es)
                j_min = ind(i_min);
                
                % Determine indices of all kinematically un-accessible 
                % solutions
                ind_isnan = find(isnan(qhat_max));
                j(ind_isnan) = false;
                
                if any(ind_isnan>j_min)
                    k = find(ind_isnan>j_min);
                    j(ind_isnan(k(1)):end) = false;
                end

                if any(ind_isnan<j_min)
                    k = find(ind_isnan<j_min);
                    j(1:ind_isnan(k(end))) = false;
                end
               
            else

                j = true(size(qhat_max));

            end
    
            
            c = [sys.Gamma_Scale,r_ana(j);
                sum(j),qhat_syn(j)/(sys.Gamma_Scale*sys.qref)];
        
            % Study asymptotic and practical stability of tuned system
            [qhat_practically_stable,qhat_stable,~,~,~,~,~,~,~] = ...
                 StabilityAnalysisLSR(c,sys,sol, ...
                            exc,'single','tuned',true,true);
            
            % Locally asymptotically stable solutions found?
            if any(~isnan(qhat_practically_stable))
                
                MaximumAmplitude(i) = max(qhat_practically_stable)/ ...
                    sys.qref;
                
            end
        end
    end
end

% Get optimum clearance
[~,~,Gamma_opt] = TunedBackbone(sys,'stable');

save([savepath 'MaximumAmplitude.mat'],'MaximumAmplitude')
save([savepath 'GammaScaleBif.mat'],'GammaScaleBif')
save([savepath 'kappa_c.mat'],'kappa_c')


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

figure(2);
hold on;
plot(kappa_c,MaximumAmplitude,'LineWidth',1.5,'Color',color.ies)
box on;
axis tight;
xlabel('$\kappa_\mathrm{c}$')
ylabel('$\hat{q}/ \hat{q}_\mathrm{ref}$')
set(gca,'XScale','log')
savefig([savepath 'MaximumAmplitude.fig'])
