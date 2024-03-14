function [qhat_max,qhat_max_violated,r] = LocalizedFrequencyAmplitudeCurve(c,sys,exc,disorder)
%LOCALIZEDFREQUENCYAMPLITUDECURVE Determine the different parts of
% localized frequency amplitude curve from contour estimation c
% c - contour of ESIM

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
qhat = abs(RecoverCondensedDOFs(sys,exc,r(indices),xi(indices),disorder));

% Maximum amplitudes
qhat_max = nan(1,length(r));
qhat_max(indices) = max(qhat,[],1);

% Index of points that fulfill kinematic constraint
qhat_max_check = true(1,length(r));
switch disorder
    case 'tuned'
        qhat_max_check(indices) = all(qhat(2:end,:)<...
            (sys.Gamma_Scale*sys.qref),1);
    case 'mistuned'
        qhat_max_check(indices) = all(qhat(2:end,:)<...
            sys.Gamma_mt(3:2:end),1);
    otherwise
        error('Case not defined.')
end

% Amplitudes with violated kinematic constraint
qhat_max_violated = nan(1,length(r));
qhat_max_violated(~qhat_max_check) = qhat_max(~qhat_max_check);

% Amplitudes with fulfilled kinematic constraint
qhat_max(~qhat_max_check) = NaN;


end

