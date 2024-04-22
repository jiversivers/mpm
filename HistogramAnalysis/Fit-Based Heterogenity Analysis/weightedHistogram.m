function [N, CENTERS] = weightedHistogram(values, weight, M)
%% weightedHistogram 
% This function will return weighted histogram counts and bin centers for
% input values and weights. M is the number of bins to bin the data into.

%% weightedHist
%
if ~exist('weight', 'var')
    weight = ones(size(values));
end

if ~exist('M', 'var')
    M = 100;
end


% Determine bins
mi = min(values, [], 'all');
ma = max(values, [], 'all');
ra = ma-mi;
wi = ra/M;
hw = wi/2;
CENTERS = mi+hw:wi:ma-hw;
    

% Ignore NaN values
vals=values(~isnan(values));
wats=weight(~isnan(values));

% Caclculate value in bin
for ii = CENTERS
    N(CENTERS==ii) = sum(wats(vals<=ii+hw));
    wats = wats(vals>ii+hw);
    vals = vals(vals>ii+hw);
end

end