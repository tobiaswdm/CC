%% Build the system

% Tuned
[sys,exc] = BuildSystem(sys,exc,'tuned');

Phi_hat_0 = [60*rand(1)*exp(2*pi*rand(1)*1i);
    220*exp(2*pi*rand(1)*1i);...
    60*rand(1)*exp(2*pi*rand(1)*1i)];

Omega_0 = 0.0037;


[Phi_hat,Omega] = ...
    SolveSlowFlowEquation(Phi_hat_0,Omega_0,1,sys,exc)

%slowfft = fft(QH(1,:))/length(QH(1,:));

%figure; plot(abs(slowfft))

%% Evaluate ESIM of GSAPR

% Auxilliary Variable
rho = (2/pi) * (1-sys.eN) / (1+sys.eN);
% Minimum amplitude - turning point of SIM
% Including safety of 0.1% higher amplitude
xi_min = 1.001 * rho/sqrt(1+rho^2);

% Clearance normalized amplitudes and frequencies
r = linspace(0.85,1.05,2000)*sys.r_k(exc.k+1);
xi = logspace(log10(xi_min),log10(20),3000);

% ESIM GSAPR
[Gamma_Scale,Xi,R] = AllSectorsESIM(xi,r,sys,exc);

% Get contour
xi_level = contourc(xi,r,Gamma_Scale,sys.Gamma_Scale*[1 1]);
xi_level(1,xi_level(2,:)==floor(xi_level(2,:))) = NaN;

% Linear FRFs
q_fixed = ComputeLinearResponse(r,sys,exc,'tuned','fixed_absorbers');
q_fixed = abs(q_fixed(1,:));
q_removed = ComputeLinearResponse(r,sys,exc,'tuned','removed_absorbers');
q_removed = abs(q_removed(1,:));


figure(2)
hold on;
plot(r,q_fixed/sys.qref,...
    'LineWidth',.5,'Color',color.reference,'DisplayName', ...
    'Fixed abs.')
plot(r,q_removed/sys.qref,'-.',...
    'LineWidth',.5,'Color',color.reference,'DisplayName', ...
    'Removed abs.')
plot(xi_level(2,:),xi_level(1,:)*sys.Gamma_Scale,'LineWidth',...
    1.5,'Color',color.ies,'DisplayName','GSAPR')
set(gca,'YScale','log')
scatter(exc.harmonic.r,sqrt(Phi_hat' * Phi_hat)/sys.qref,...
    50,'filled','DisplayName','AQPR')
axis tight;
xlabel('$r$')
ylabel('$\hat{q}/\hat{q}_\mathrm{ref}$')
legend;
box on;