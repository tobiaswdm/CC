function N_SIPP = CountImpacts(UA,N_T,mode)

switch mode
    case 'average'
        N_SIPP = zeros(size(UA,1),1);
        
        for i = 1:(size(UA,2)-1)
            
            simp = sign(UA(:,i)) ~= sign(UA(:,i+1)); % Find indices of absorbers with significant impact
            
            N_SIPP(simp) = N_SIPP(simp) + 1;
            
        end
        
        N_SIPP = N_SIPP/N_T;
    
    case 'convolution'
        % Count number of impacts (incuding insignificant) in each period
        
        % Number of points per period
        N_P = (size(UA,2)-1)/N_T;

        N_SIPP = zeros(size(UA,1),N_T);

        for j = 0:(N_T-1)

            for i = 1:N_P
                
                % Detect change in absorber velocity -> impact
                imp = UA(:,j*N_P+i) ~= UA(:,j*N_P+i+1);
                
                % Check if impact is not accurately resolved and
                % takes more than one time step to avoid counting it
                % twice
                if i>1
                    imp = imp & (UA(:,j*N_P+i-1) == UA(:,j*N_P+i));
                end
                
                N_SIPP(imp,j+1) = N_SIPP(imp,j+1) + 1;

            end

        end

    otherwise
        error('Case not defined.')
end

end