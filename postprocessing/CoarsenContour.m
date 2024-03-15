function [c_coarse] = CoarsenContour(c,stepsize)
%COARSENCONTOUR Coarsen the low level contour estimation to include in
% numerical stability analysis
%
% c - low level contour for single clearance
% stepsize - consider only every stepsize-th point

% Capture separation points between isolated branches
indices = floor(c(2,:))==c(2,:);

% Set other points in coarser sampling to true
indices(2:stepsize:length(indices)) = true;

% Coarse contour estimation
c_coarse = c(:,indices);

end

