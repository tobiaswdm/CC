function y = ThreeParamWeibull(x,a,b,c,type)
%THREEPARAMWEIBULL Summary of this function goes here
%   Detailed explanation goes here


switch type
    case 'pdf' % PDF
        x = x-c;
        x(x<0) = 0;
        y = wblpdf(x,a,b);
    case 'cdf' % CDF
        x = x-c;
        x(x<0) = 0;
        y = wblcdf(x,a,b);
    case 'log' % Log PDF
        x = x-c;
        x(x<0) = 0;
        y = log(b/a) + (b-1)*log(x/a) - (x/a).^b;
    case 'nlog' % Negative Log Likelyhood
        x = x-c;
        x(x<0) = 0;
        y = -log(b/a) - (b-1)*log(x/a) + (x/a).^b;
    case 'inv'
        y = wblinv(x,a,b);
        y = y + c;
    otherwise
        error('Case not defined.')
end

end

