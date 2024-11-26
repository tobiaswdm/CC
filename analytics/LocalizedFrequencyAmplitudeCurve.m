%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CC (pronounced Sisi) is a tool that performs numerical and/or
% analytical analyses on a Cylic Chain of Oscialltors with Vibro-Impact
% Nonlinear Energy Sinks (VI-NESs)
%
% The Code for CC was written by:
% Tobias Weidemann - (C) 2024
% University of Stuttgart, Germany
% Institute of Aircraft Propulsion Systems
%
% Contact: tobias.weidemann@ila.uni-stuttgart.de
%
% Feel free to use, share and modify under the GPL-3.0 license.
% CC is purely academic and comes with no warranty.
% If you use CC for your own research, please refer to the paper:
%
% T. Weidemann, L. A. Bergman, A. F. Vakakis, M. Krack. (2024)
% "Energy Transfer and Localization in a Forced Cyclic Chain of
% Oscillators with Vibro-Impact Nonlinear Energy Sinks".
% Manuscript submitted to Nonlinear Dynamics
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [qhat_max,qhat_max_violated,qhat_localized,r] = ...
    LocalizedFrequencyAmplitudeCurve(c,sys,exc,pattern,disorder)
%LOCALIZEDFREQUENCYAMPLITUDECURVE Determine the different parts of
% localized frequency amplitude curve from contour estimation c
% c - contour of FRS
% pattern 'single' or 'opposing'

% Find beginning of level curves in c
c(2,floor(c(2,:))==c(2,:)) = NaN;

% Frequencies
r = c(1,2:end);
% Ampltiudes in reference sector
xi = c(2,2:end);

% Determine indices of actual level curves
% without NaN spacers
indices = ~isnan(xi);

% Recover real-valued amplitudes of all sectors
[qhat,~,~] = RecoverCondensedDOFs(sys,exc,r(indices),xi(indices),pattern, ...
    disorder);
qhat = abs(qhat);

% Maximum amplitudes
qhat_max = nan(1,length(r));
qhat_max(indices) = max(qhat,[],1);

% Which sectors are non-synchronized
switch pattern
    case 'single'
        non_synchronized = 2:sys.N_s;
    case 'opposing'
        non_synchronized = 1:sys.N_s;
        non_synchronized([1,sys.N_s/2+1]) = [];
end

% Index of points that fulfill kinematic constraint
qhat_max_check = true(1,length(r));
switch disorder
    case 'tuned'
        qhat_max_check(indices) = all(qhat(non_synchronized,:)<...
            (sys.Gamma_Scale*sys.qref),1);
        qhat_localized = sys.Gamma(1)*xi;
    case 'mistuned'
        qhat_max_check(indices) = all(qhat(non_synchronized,:)<...
            sys.Gamma_mt(3:2:end),1);
        qhat_localized = sys.Gamma_mt(1)*xi;
    otherwise
        error('Case not defined.')
end

% Amplitudes with violated kinematic constraint
qhat_max_violated = nan(1,length(r));
qhat_max_violated(~qhat_max_check) = qhat_max(~qhat_max_check);

% Amplitudes with fulfilled kinematic constraint
qhat_max(~qhat_max_check) = NaN;

end

