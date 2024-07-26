function sol = ConfigureIntegrator( ...
    sol,sys,exc,initial_condition,cut_transient,disorder)
% Configure the solver parameters depending on the system and task
%
% initial_condition - 'zero', 'random' or 'localized'
%                      if any other string is used than the ICs won't be
%                      changed from the input struct sol
% cut_transient - true or false (simulates through transient parst first
%                                before outputting)
% disorder - 'tuned' for tuned system and anything else for disordered
%                   system
    
    % Assign initial condition
    switch initial_condition
        case 'zero'
            sol.q0 = zeros(sys.N_s,1);
            sol.u0 = zeros(sys.N_s,1);
            sol.ua0 = zeros(sys.N_s,1);
            sol.qa0 = zeros(sys.N_s,1);
            % Place NES at cavity walls
            if strcmp(disorder,'tuned')
                sol.qa0 = -0.99*sys.Gamma(1:2:2*sys.N_s)+sol.q0;
            else
                sol.qa0 = -0.99*sys.Gamma_mt(1:2:2*sys.N_s)+sol.q0;
            end
        case 'random'
            % Chose random to break symmetry of the system
            q0 = 1e-4 * sys.qref * (2*rand(sys.N_s,1)-1);
            % Moreau solver uses modal coordinates for 
            sol.q0 = transpose(sys.Phi)*q0;
            sol.u0 = zeros(sys.N_s,1);
            sol.ua0 = zeros(sys.N_s,1);
            if strcmp(disorder,'tuned')
                sol.qa0 = -0.99*sys.Gamma(1:2:2*sys.N_s)+q0;
            else
                sol.qa0 = -0.99*sys.Gamma_mt(1:2:2*sys.N_s)+q0;
            end
        case 'localized'
            % Chose localized to assign one sector with
            % a large initial displacement
            q0 = zeros(sys.N_s,1);
            q0(1) = 2*abs(sys.qref);
            sol.q0 = transpose(sys.Phi)*q0;
            sol.u0 = zeros(sys.N_s,1);
            sol.ua0 = zeros(sys.N_s,1);
            if strcmp(disorder,'tuned')
                sol.qa0 = -0.99*sys.Gamma(1:2:2*sys.N_s)+q0;
            else
                sol.qa0 = -0.99*sys.Gamma_mt(1:2:2*sys.N_s)+q0;
            end
    end

    % Set up time domain
    Tau = 2*pi/exc.harmonic.r;  % Period length
    if cut_transient
        % Determine when linear transient response would be decayed
        % by 99.9 percent
        tau_decay = log(1000)/sys.D;
        % Corresponsed to how many excitation periods?
        sol.NP_trans = sol.N_P*ceil(tau_decay/Tau);
    else
        sol.NP_trans = 0;
    end
    
    sol.dtau = Tau/sol.N_P;
    sol.N_Save = sol.N_P/sol.N_Sample;
    sol.N_tau = sol.N_Tau*sol.N_P;

end
