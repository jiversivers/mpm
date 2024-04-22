function [curveVals, fitStats, R2] = freqCurve(CENTERS, N)

%% freqCurve
% This function fits the input histogram with a 1-5 component Gaussian
% Mixture Model (GMM) The Akaike Information Criterion local minimum is
% used to determine the best fit for the data from each fit tested (1-5
% Gaussians).

% Updated: 11-Aug-2022 by Jesse Ivers
% Recoded histogram function and broke off of this function to clean up.

%% FITTING
[xData, yData] = prepareCurveData(CENTERS, N); 
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';

% Try 1-5 Components
gmmFit = cell(1,5);
gof = cell(1,5);
RSS = [sum(N.^2), zeros(1,5)];

for popN = 1:5
    ft = fittype(['gauss', num2str(popN)]);
    opts.Lower = [opts.Lower, 0 0 0];
    opts.Upper = [opts.Upper, Inf 1 1];

    % In case of failed fit at given component count
    try
        [gmmFit{popN}, gof{popN}] = fit(xData, yData, ft, opts);
        RSS(popN+1) = gof{popN}.sse;
    catch
        gmmFit{popN} = missing;
        gof{popN} = missing;
        RSS(popN+1) = empty;
    end
end

% Calculate fit scores and choose
n = numel(CENTERS);
p = 3*(0:5);
AIC = 2*p + n*(log(RSS/n)) + (2*p.*(p+1))./(n-p-1);
dAIC = AIC(2:end)-AIC(1:end-1);
ddAIC = dAIC(2:end)-dAIC(1:end-1);
% totalCount = find(AIC == min(AIC))-1;     % Min AIC
 totalCount = find(ddAIC > 0, 1)+1;          % First inflection point
% totalCount = find(dAIC > 0, 1);  % Greatest decrease in rate of improvement
 if isempty(totalCount)
     totalCount = find(dAIC == max(dAIC));
 end


% Get individual curve parameters
[mu, sig, port, curveVals] = gaussStats(gmmFit{totalCount}, CENTERS);

% Sort final outputs by ascending means
[mu, I] = sort(mu);
sig = sig(I);
port = port(I);
curveVals = curveVals(I,:);

% Final outputs
fitStats = [mu, sig, port];
R2 = gof{totalCount}.rsquare;

end