function [Q,QA,U,UA,T] = MoreauIntegration(sys,exc,sol,disorder)

% Assign Variables
Ma = sys.Ma;
eN = sys.eN;
fhat = exc.F_mod;
WN = sys.WN_mod;
WNA = sys.WNA;

q0 = sol.q0;
qa0 = sol.qa0;
u0 = sol.u0;
ua0 = sol.ua0;

N_transient = sol.NP_trans;
dt = sol.dtau;
N_Save = sol.N_Save;
Nt = sol.N_tau;

switch disorder
    case 'tuned'

        M = sys.mu;
        C = sys.beta;
        K = sys.kappa;
        g = sys.Gamma;

    case 'mistuned'

        M = sys.mu_mt;
        C = sys.beta_mt;
        K = sys.kappa_mt;
        g = sys.Gamma_mt;

    otherwise
        error('Case not defined.')
end


% dt - Timestep
% N_transient - Time steps that are ignored
% Nt - Time steps of stationary response
% N_Save - Save every N_Save-th time step of stationary response

if rem(Nt,N_Save)
    error('Saving time steps must be multiple of simulation time steps.')
end

L = Nt/N_Save+1;

% Set up matrices
n = sys.N_s;            % Number of DOF in linear system (2*n states)
nA = length(qa0);       % Number of absorbers
Q = zeros(n,L);         % Initialize Displacement
QA = zeros(nA,L);       % Initialize Displacement Absorber
U = zeros(n,L);         % Initialize Velocity
UA = zeros(nA,L);       % Initialize Velocity Absorber
T = zeros(1,L);         % Time

% Set inital values
QB = q0;
QBA = qa0;
UB = u0;
UBA = ua0;
t = 0;
saveindex = N_Save;
k = 2;
        
% Symmetric moreau as derived in
% 'Comparison of Moreau-type integrators based on the time finite
% element discretization of the virtual action' - Capobianco et al., ENOC 2017
% This paper uses a different (confusing) indexing
% k-1 = beginning of timestep
% k = midpoint
% k+1 = end of timestep

% Set up matrices
MC = M + 0.5*dt*C;       
MC_inv = MC\eye(n);
Ma_inv = Ma\eye(nA);
MC_invM = MC_inv*M;
for j = 1:N_transient
    
    % Approximate Q at midpoint
    QM = QB + 0.5*dt*UB; % Masses
    QMA = QBA + 0.5*dt*UBA; % Absorbers
    tM = t + 0.5*dt;
    
    % Collective symmetric forces at midpoint
    hM = ExcitationTime(tM,fhat,exc)-K*QM-0.5*C*UB;
    
    % Contact detection at midpoint (WN,WNA = const.)
    gN = g + WN'*QM + WNA'*QMA;
    IC = find(gN<=0); nC = length(IC);
    
    % Handle contacts
    if nC>0
        % Set up inclusion problem
        
        % Normalized force directions of active contacts
        WN_active = WN(:,IC); WNA_active = WNA(:,IC);
        % Contact velocity at beginning of timestep
        
        % Initialize contact efforts
        P_prev = zeros(nC,1); % Cold start (all zero)

        
        iter = 1;
        converged = false;
        
        % Delassu Matrix
        G = WN_active' * MC_inv * WN_active + ...
            WNA_active' * Ma_inv * WNA_active;
        r = sol.alpha./diag(G);
        
        rWUB = eN*r.*(WN_active'*UB + WNA_active'*UBA);
        UE = MC_invM*UB + MC_inv*(hM*dt+WN_active*P_prev);
        UEA = UBA + Ma_inv*(WNA_active*P_prev);
        
        while ~converged
            switch sol.solver
                case 'JOR'
                    P = -min(-P_prev + rWUB + r.*(WN_active'*UE + WNA_active' * UEA),0);
                    UE = UE + MC_inv*WN_active*(P-P_prev);
                    UEA = UEA + Ma_inv*WNA_active*(P-P_prev);
                case 'fsolve'                            
                    % Set up function, where x = [UE;UEA;PP]
                    % UE = x(1:n); UEA = x((n+1):(n+nA)); PP = x((n+nA+1))
                    nUE_start = n+1;
                    nUE_end = n+nA;
                    nPP_start = n+nA+1;
                    nPP_end = n+nA+nC;
                    f = @(x) [MC*x(1:n)-M*UB - hM*dt - WN_active*x(nPP_start:nPP_end);...
                        Ma*(x(nUE_start:nUE_end)-UBA) - WNA_active*x(nPP_start:nPP_end);...
                        x(nPP_start:nPP_end) + min(-x(nPP_start:nPP_end) + rWUB +...
                        sol.rN*(WN_active'*x(1:n) + WNA_active' * x(nUE_start:nUE_end)),0)];
                    options = optimoptions('fsolve','MaxIterations',sol.maxiter,'Display','off');
                    [xout,~,exitflag,~] = fsolve(f,[UB;UBA;P_prev],options);
                    UE = xout(1:n);
                    UEA = xout(nUE_start:nUE_end);
                    P = xout(nPP_start:nPP_end);
                    if exitflag == 0 || exitflag == -2
                        warning('No convergence in time step.')
                    end
                    P_prev = P;                               
                otherwise
                    error(['Unknown solver ' sol.solver '.']);
            end
            err = sum(abs(P-P_prev));
            P_prev = P;
            converged = err <= sol.tol;
            
            iter = iter+1;
            if iter > sol.maxiter
                converged = true;
                warning('No convergence in time step.')
            elseif any(isnan([UE;UEA]))
                error('NaN in prox iteration.')
            end
        end
        
    else
        UE = MC_invM*UB + MC_inv*hM*dt;
        UEA = UBA;
    end
    
    t = t+dt;
    QE = QM + 0.5*dt*UE;
    QEA = QMA + 0.5*dt*UEA;
    QB = QE;
    QBA = QEA;
    UB = UE;
    UBA = UEA;           
end

Q(:,1) = QB;      % Initialize Displacement
QA(:,1)= QBA;     % Initialize Displacement Absorber
U(:,1) = UB;      % Initialize Velocity
UA(:,1) = UBA;    % Initialize Velocity Absorber
T(1) = t;         % Time
for j = 1:Nt
    
    % Approximate Q at midpoint
    QM = QB + 0.5*dt*UB; % Masses
    QMA = QBA + 0.5*dt*UBA; % Absorbers
    tM = t + 0.5*dt;
    
    % Collective symmetric forces at midpoint
    hM = ExcitationTime(tM,fhat,exc)-K*QM-0.5*C*UB;
    
    % Contact detection at midpoint (WN,WNA = const.)
    gN = g + WN'*QM + WNA'*QMA;
    IC = find(gN<=0); nC = length(IC);
    
    % Handle contacts
    if nC>0
        % Set up inclusion problem
        
        % Normalized force directions of active contacts
        WN_active = WN(:,IC); WNA_active = WNA(:,IC);
        % Contact velocity at beginning of timestep
        
        % Initialize contact efforts
        P_prev = zeros(nC,1); % Cold start (all zero)

        
        iter = 1;
        converged = false;
        
        % Delassu Matrix
        G = WN_active' * MC_inv * WN_active + ...
            WNA_active' * Ma_inv * WNA_active;
        r = sol.alpha./diag(G);
        
        rWUB = eN*r.*(WN_active'*UB + WNA_active'*UBA);
        UE = MC_invM*UB + MC_inv*(hM*dt+WN_active*P_prev);
        UEA = UBA + Ma_inv*(WNA_active*P_prev);
        
        while ~converged
            switch sol.solver
                case 'JOR'
                    P = -min(-P_prev + rWUB + r.*(WN_active'*UE + WNA_active' * UEA),0);
                    UE = UE + MC_inv*WN_active*(P-P_prev);
                    UEA = UEA + Ma_inv*WNA_active*(P-P_prev);
                case 'fsolve'                            
                    % Set up function, where x = [UE;UEA;PP]
                    % UE = x(1:n); UEA = x((n+1):(n+nA)); PP = x((n+nA+1))
                    nUE_start = n+1;
                    nUE_end = n+nA;
                    nPP_start = n+nA+1;
                    nPP_end = n+nA+nC;
                    f = @(x) [MC*x(1:n)-M*UB - hM*dt - WN_active*x(nPP_start:nPP_end);...
                        Ma*(x(nUE_start:nUE_end)-UBA) - WNA_active*x(nPP_start:nPP_end);...
                        x(nPP_start:nPP_end) + min(-x(nPP_start:nPP_end) + rWUB +...
                        sol.rN*(WN_active'*x(1:n) + WNA_active' * x(nUE_start:nUE_end)),0)];
                    options = optimoptions('fsolve','MaxIterations',sol.maxiter,'Display','off');
                    [xout,~,exitflag,~] = fsolve(f,[UB;UBA;P_prev],options);
                    UE = xout(1:n);
                    UEA = xout(nUE_start:nUE_end);
                    P = xout(nPP_start:nPP_end);
                    if exitflag == 0 || exitflag == -2
                        warning('No convergence in time step.')
                    end
                    P_prev = P;                               
                otherwise
                    error(['Unknown solver ' sol.solver '.']);
            end
            err = sum(abs(P-P_prev));
            P_prev = P;
            converged = err <= sol.tol;
            
            iter = iter+1;
            if iter > sol.maxiter
                converged = true;
                warning('No convergence in time step.')
            elseif any(isnan([UE;UEA]))
                error('NaN in prox iteration.')
            end
        end
        
        
    else
        UE = MC_invM*UB + MC_inv*hM*dt;
        UEA = UBA;
    end
    
    t = t+dt;
    QE = QM + 0.5*dt*UE;
    QEA = QMA + 0.5*dt*UEA;
    QB = QE;
    QBA = QEA;
    UB = UE;
    UBA = UEA;
    
    if j ==  saveindex
        Q(:,k) = QE;      % Initialize Displacement
        QA(:,k)= QEA;     % Initialize Displacement Absorber
        U(:,k) = UE;      % Initialize Velocity
        UA(:,k) = UEA;    % Initialize Velocity Absorber
        T(k) = t;         % Time
        k = k+1;
        saveindex = saveindex + N_Save;
    end
end
               
end



function [f_out] = ExcitationTime(t,F,exc)
% Calculate excitation force at given time t

% t - time instance or 1xNt vector
% F - complex forcing vector
% exc - struct with excitation data

switch exc.type
    case 'harmonic' % Harmonic excitation
        f_out = real(F*exp(1i*exc.harmonic.r*t));
    case 'sweep'
        f_out = real(F*exp(1i*(exc.sweep.r0 + ...
            0.5*(exc.sweep.re-exc.sweep.r0)/exc.sweep.tau*t)*t));
    case 'transient'
        f_out = 0*F;
    otherwise
        error('Excitation not defined.')
end

end
