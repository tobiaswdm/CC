function [E] = ModalEnegies(sys,sol,ETA,CHI)
%MODALENEGIES Summary of this function goes here
%   Detailed explanation goes here

movRMS = dsp.MovingRMS(sol.N_Sample);

ETAsq = (movRMS(ETA')').^2;
CHIsq = (movRMS(CHI')').^2;

E = zeros(sys.k(end)+1,size(ETA,2));

E(1,:) = 0.5*(CHIsq(1,:)+sys.r_k(1)^2 * ETAsq(1,:));

for k = 1:(floor(sys.N_s/2)-1)

    E(k+1,:) = 0.5*(CHIsq(2*k,:)+sys.r_k(1)^2 * ETAsq(2*k,:) + ...
                CHIsq(2*k+1,:)+sys.r_k(1)^2 * ETAsq(2*k+1,:));

end

if rem(sys.N_s,2)
    k = floor(sys.N_s/2);
    E(k+1,:) = 0.5*(CHIsq(2*k,:)+sys.r_k(1)^2 * ETAsq(2*k,:) + ...
                CHIsq(2*k+1,:)+sys.r_k(1)^2 * ETAsq(2*k+1,:));
else
    E(end,:) = 0.5*(CHIsq(sys.N_s,:)+sys.r_k(1)^2 * ETAsq(sys.N_s,:));
end

