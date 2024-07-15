function [E,E_avg] = ModalEnegies(sys,sol,ETA,CHI)
%MODALENEGIES Determine the instantaneus modal energies E and their average
% over time E_avg.

movRMS = dsp.MovingRMS(sol.N_Sample);

ETAsq = (movRMS(ETA')').^2;
CHIsq = (movRMS(CHI')').^2;

E = zeros(sys.k(end)+1,size(ETA,2));
E_avg = zeros(sys.k(end)+1,1);

E(1,:) = 0.5*((1-sys.epsilon_a)*CHIsq(1,:)+sys.r_k(1)^2 * ETAsq(1,:));
E_avg(1) = mean(0.5*((1-sys.epsilon_a)*CHI(1,:).^2+...
    (sys.r_k(1) .* ETA(1,:)).^2));

for k = 1:(floor(sys.N_s/2)-1)

    E(k+1,:) = 0.5*((1-sys.epsilon_a)*(CHIsq(2*k,:)+CHIsq(2*k+1,:))+...
        sys.r_k(k+1)^2 * (ETAsq(2*k,:) + ETAsq(2*k+1,:)));

    E_avg(k+1) = mean(0.5*((1-sys.epsilon_a)*(CHI(2*k,:).^2+...
        CHI(2*k+1,:).^2)+sys.r_k(k+1)^2 * (ETA(2*k,:).^2 +...
        ETA(2*k+1,:).^2)));

end

if rem(sys.N_s,2)
    k = floor(sys.N_s/2);
    E(k+1,:) = 0.5*((1-sys.epsilon_a)*(CHIsq(2*k,:)+CHIsq(2*k+1,:))+...
        sys.r_k(k+1)^2 * (ETAsq(2*k,:) + ETAsq(2*k+1,:)));
    E_avg(end) = mean(0.5*((1-sys.epsilon_a)*(CHI(2*k,:).^2+CHI(2*k+1,:).^2)+...
        sys.r_k(k+1)^2 * (ETA(2*k,:).^2 + ETA(2*k+1,:).^2)));
else
    E(end,:) = 0.5*((1-sys.epsilon_a)*CHIsq(sys.N_s,:)+...
        sys.r_k(end)^2 * ETAsq(sys.N_s,:));
    E_avg(end) = mean(0.5*((1-sys.epsilon_a)*CHI(sys.N_s,:).^2 +...
        sys.r_k(end)^2 * ETA(sys.N_s,:).^2));
end

end

