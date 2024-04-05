function ProjectOnSIM(Q,QA,sol,sys,color,savepath,disorder)
%PROJECTONSIM Porject motion of first sector onto stable branch of 
% SIM derived by Gendelman
%
% Q - [1,N_Tau * N_Sample] Displacement Host Oscillator
% QA - [1,N_Tau * N_Sample] Displacement Absorber

% Check if length is consistent
if length(Q) ~= sol.N_Sample*sol.N_Tau
    error(['Length of array must be integer multiple of' ...
        ' excitation frequency.'])
end

% Time Vector
Periods = 1:sol.N_Tau;
% Oscillator amplitude
qhat = zeros(1,sol.N_Tau);
% Absorber amplitude
qahat = zeros(1,sol.N_Tau);

% Extract amplitude values in each excitaiton period
for i = 1:sol.N_Tau
    
    % Oscillator amplitude
    qhat(i)=0.5*(max(Q(((i-1)*sol.N_Sample+1):(i*sol.N_Sample)),[],2)-...
        min(Q(((i-1)*sol.N_Sample+1):(i*sol.N_Sample)),[],2));
    
    % Absorber amplitude
    qahat(i)=0.5*(max(QA(((i-1)*sol.N_Sample+1):(i*sol.N_Sample)),[],2)-...
        min(QA(((i-1)*sol.N_Sample+1):(i*sol.N_Sample)),[],2));

end

% Normalize by clearance
switch disorder
    case 'tuned'
        xi = qhat/sys.Gamma(1);
        theta = qahat/sys.Gamma(1);
    case 'mistuned'
        xi = qhat/sys.Gamma_mt(1);
        theta = qahat/sys.Gamma_mt(1);
    otherwise
        error('Case not defined.')
end


%% Determine analytical branch

% Auxilliary Variable
rho = (2/pi) * (1-sys.eN) / (1+sys.eN);

% Minimum norm. amplitude
xi_min = 1.0001 * rho/sqrt(1+rho^2);
% Vector for plotting analytical sim
xi_ana = linspace(xi_min,1.2*max(xi),2000);
% Grid for surface plot of SIM
[TAU_ANA,XI_ANA] = meshgrid([1,sol.N_Tau],xi_ana);
% Absorber Amplitude
THETA_ANA = repmat((1+sqrt((1+rho^2)*(xi_ana').^2 - rho^2))/...
    (1+rho^2),[1,2]);

% Plot motion on SIM
figure;
surf(TAU_ANA,XI_ANA,THETA_ANA,'EdgeAlpha',0,'FaceColor',color.show,...
    'FaceAlpha',0.2)
hold on;
plot3(Periods,xi,theta,'LineWidth',2,'Color',color.analytics)
plot3(Periods,xi,THETA_ANA(1,1)*ones(1,sol.N_Tau),...
    'LineWidth',2,'Color',color.background)
plot3(Periods,xi_ana(end)*ones(1,sol.N_Tau),theta,...
    'LineWidth',2,'Color',color.background)
xlabel('$r\tau/(2\pi)$')
ylabel('$\xi$')
zlabel('$\vartheta$')
axis tight
title('Motion on SIM')
box on;
if ~isempty(savepath)
    savefig([savepath 'SIM_motion.fig'])
end



end

