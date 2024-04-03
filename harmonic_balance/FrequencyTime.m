function [transformation] = FrequencyTime(domain,resolution,method)
%FREQUENCYTIME Get Fourier coefficients for time input or time domain
% from Fourier Coefficients
%
% H - Harmonic order
%
% Freq_to_time:
% domain - [2H+1,1] Complex Fourier coefficients
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
        
        % Harmonic order
        H = (size(domain,1)-1)/2;
        
        % Check if resolution is sufficient
        if resolution < length(domain)
            error('Chose higher sampling rate.')
        end

        % harmonics
        h = -H : H;
        
        % Complex Harmonic Basis functions
        W = exp(1i*(2*pi)*(0:(resolution-1))'*h /resolution);
        
        % Transform back to timedomain
        transformation = (W*domain);

    case 'Time_to_Freq'

        % length of signal
        % assumed to be even
        L = length(domain);

        if L < (2*resolution+1)
            error('Chose higher sampling rate.')
        end

        % DFT using FFT
        freqdomain = fft(domain,[],1);

        % Get only symmetric half and select number of harmonics
        transformation = [freqdomain((end-resolution+1):end,:);...
                          freqdomain(1:(resolution+1),:)]/L;

    otherwise

        error('Case not defined.')
end


end

