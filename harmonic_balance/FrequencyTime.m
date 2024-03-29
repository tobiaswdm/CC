function [transformation] = FrequencyTime(domain,resolution,method)
%FREQUENCYTIME Get Fourier coefficients for time input or time domain
% from Fourier Coefficients
%
% H - Harmonic order
%
% Freq_to_time:
% domain - [H+1,1] Complex Fourier coefficients
% resolution - number of sampling points per period
%
% Desired resolution must be larger or equal to 2H+1
%
% Time_to_freq:
% domain - [N_time,1] Real-valued time series
% resolution - harmonic order H
%
% N_time must be larger or equal to 2H+1

switch method
    case 'Freq_to_Time'
        
        % Check if resolution is sufficient
        if resolution < (2*length(domain)-1)
            error('Chose higher sampling rate.')
        end
        
        % Complex Harmonic Basis functions
        W = exp(1i*(2*pi)*(0:(resolution-1))' * ...
            (0:(length(domain)-1))/resolution);
        
        % Transform back to timedomain
        transformation = real(W*domain);

    case 'Time_to_Freq'

        % length of signal
        % assumed to be even
        L = length(domain);

        if L < (2*resolution+1)
            error('Chose higher sampling rate.')
        end

        % DFT using FFT
        freqdomain = fft(domain,1);

        % Get only symmetric half and select number of harmonics
        transformation = [freqdomain(1);2*freqdomain(2:(resolution+1))]/L;

    otherwise

        error('Case not defined.')
end


end

