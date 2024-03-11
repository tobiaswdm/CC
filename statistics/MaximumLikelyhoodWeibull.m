function [a,b,c] = MaximumLikelyhoodWeibull(x,type)
% Maximum Likelyhood Estimation of three Parameter Weibull to data x

% Get minimum in x for bounds and start of parameter c
x_min = min(x,[],'all');

switch type
    case 'pdf'
        pdffunc = @(x,a,b,c) ThreeParamWeibull(x,a,b,c,'pdf');
        params = mle(x,'pdf',pdffunc,'Start',[2 3 0.8*x_min],...
            'Options',statset('FunValCheck','off'),...
            'LowerBound',[0 0 -Inf],'UpperBound',[Inf Inf x_min]);
    case 'nlog'
        nloglffunc =  @(x,a,b,c) ThreeParamWeibull(x,a,b,c,'nlog');
        params = mle(x,'nloglf',nloglffunc,'Start',[2 3 0.8*x_min],...
            'LowerBound',[0 0 -Inf],'UpperBound',[Inf Inf x_min]);
    case 'log'
        loglffunc =  @(x,a,b,c) ThreeParamWeibull(x,a,b,c,'log');
        params = mle(x,'logpdf',loglffunc,'Start',[2 3 0.8*x_min],...
            'LowerBound',[0 0 -Inf],'UpperBound',[Inf Inf x_min]);
    otherwise
        error('Case not defined.')
end

a = params(1);
b = params(2);
c = params(3);
end

