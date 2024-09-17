function [IPR,LF] = LocalizationMeasures(qhat,sys)

    % Compute the Localization Measures
    % 'IPR' Inverse participation ratio as defined in
    %                 https://doi.org/10.1103/PhysRevB.103.024106
    %  or 'LF' Localization Factor as defined in
    %                 https://doi.org/10.1115/1.2985074 
    %
    % Input:
    % qhat - [Ns,1] Time-averaged amplitudes
        
    % IPR
    IPR = sum(qhat.^4)/(sum(qhat.^2))^2;
    
    % LF
    chi = max(qhat)/rms(qhat);
    LF = (chi-1)/(sqrt(sys.N_s)-1);

end