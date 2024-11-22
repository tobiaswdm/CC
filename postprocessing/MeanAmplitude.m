function [AMEAN,ASTD] = MeanAmplitude(x,NP)
    % NP - Sample Points per excitation period
    % x - Reponse

    % AMEAN - mean amplitude of excitation period
    % SPEC - complex FFT of the signal
        
    % Take maximum value
    NR = size(x,2)/NP;
    n = size(x,1);
    
    if NR ~= floor(NR)
        error('Length of array must be integer multiple of excitation frequency.')
    end
    
    amplitudes = zeros(n,NR);
    
    for i = 1:NR
        %X = abs(fft(x(:,((i-1)*NP+1):(i*NP)),[],2)/NP);
        %amplitudes(:,i) = 2*X(:,2);
        amplitudes(:,i)=0.5*(max(x(:,((i-1)*NP+1):(i*NP)),[],2)-min(x(:,((i-1)*NP+1):(i*NP)),[],2));
    end
    
    AMEAN = mean(amplitudes,2);
    ASTD = std(amplitudes,1,2);

end