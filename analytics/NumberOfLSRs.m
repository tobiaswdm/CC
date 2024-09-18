function [N_LSR_tuned,N_LSR_mistuned] = NumberOfLSRs(N_s)
%NUMBEROFLSRS Returns the number of possible coexisting unique LSRs
%   N_s - Vector of Number of sectors
%
%   N_LSR_tuned - Number of unqiue coexisiting LSRs in tuned system
%   N_LSR_mistuned - Number of unqiue coexisiting LSRs in mistuned system

% Each sector is either synchronized or non-synchronized
% Two substracted to exclude synchronization in all (GSR) or no sectors
% (linear response)
N_LSR_mistuned = 2.^N_s - 2;

% Each sector is either synchronized or non-synchronized
% Two substracted to exclude synchronization in all (GSR) or no sectors
% (linear response)
% Necklace counting problem: 
% https://en.wikipedia.org/wiki/Necklace_(combinatorics)

N_LSR_tuned = zeros(size(N_s));
for i = 1:numel(N_s)
    N_LSR_tuned(i)  = sum(2.^gcd(1:N_s(i),N_s(i)))/N_s(i) - 2;
end

end

