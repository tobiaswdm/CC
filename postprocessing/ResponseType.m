function type = ResponseType(nsipp)

% Determine the type of response based on the number of significant
% Impacts per excitation period in the system
% 
% type = 0 - GSR
% type = 1 - LSR
% type = 2 - SMR (No differentiation between SSMR and ASMR!)

type = nan(1,size(nsipp,2));

% GSR
index = min(nsipp,[],1)>=1.98; % All sectors in 1:1 resonance
type(index) = 0;

% LSR
% At least 1 sector in 1:1 resonance and at least 1 sector
% not in 1:1 resonance
index = max(nsipp,[],1)>=1.98 & ...
        min(nsipp,[],1)<1.98;
type(index) = 1;

% SMR
% No sector in 1:1 rsonance but impacts in all sectors
index = max(nsipp,[],1)<1.98 & ...
        min(nsipp,[],1)>0;
type(index) = 2;

if any(isnan(type))
    warning('Response type could not be indentified.')
end

end