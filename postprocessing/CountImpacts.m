function N_SIPP = CountImpacts(UA,N_T)

N_SIPP = zeros(size(UA,1),1);

for i = 1:(size(UA,2)-1)
    
    simp = sign(UA(:,i)) ~= sign(UA(:,i+1)); % Find indices of absorbers with significant impact
    
    N_SIPP(simp) = N_SIPP(simp) + 1;
    
end

N_SIPP = N_SIPP/N_T;

end