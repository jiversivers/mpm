function [fig, subMean, subVar, z, p] = WeightedHistFits(val, int, pltTitle, xaxTitle)

%% weightedHistFits
% This function creates an intensity-weighted histogram of values in an
% image. This histogram is then fit with Gaussian curve(s). The raw
% histogram and the fit are both plotted automattically. The Akaike
% Information Criterion local minimum is used to determine the best fit for
% the data from 1 to 8 Gaussians. This fit is then deconstructed into its
% individual Gaussian curves which are then compared through a simplified
% z-test. This z-test takes the square difference of the means divided by
% the square root of the sum of the variances. From this z score, a p value
% is determined. For those pairs of individual curves with a p value of no
% more than 0.05, the associated individual curves are added to the plot,
% along with markers for their corresponding means and FWHM values.
% Alongside this plot, mean and variance for all Gaussian curves in the
% model are listed, along with an annotation of the number of components in
% the fit.

% Updated: 26-Jul-2022 by Jesse Ivers

% Calculate weighted histogram (column-wise cat, so values are still
% aligned)
histVals = val(~isnan(val));
weights = int(~isnan(val));
[N, binEdges, binInd] = histcounts(histVals,50);
histOut = accumarray(binInd, weights);

% Find center of each bin
for ii = 1:length(N)
    binCenter(ii)= (binEdges(ii) + binEdges(ii+1))/2;
end

% Ignore dark and sat pixel bins and normalize weight to sum to 1
binCenter = binCenter(2:end-1);
histOut = histOut(2:end-1)./sum(histOut(2:end-1),'all');

%% Fitting
[xData, yData] = prepareCurveData(binCenter, histOut');

% Set up fit options
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';

% Fit data with 1-8 Gaussians
rsquare = [0];
fitresult = cell(8,1);

for popN = 1:8
    try
        ft = ['gauss', num2str(popN)];
        [fitresult{popN}, gof] = fit(xData, yData, ft, opts);
        rsquare(popN) = gof.rsquare;
        RSS(popN) = gof.sse;
    catch
        rsquare(popN) = 0;
        RSS(popN) = inf;
    end
end

% Determine best fit
% Calculate AIC Score for each fit
T = 48;                                                 % Sample size of histogram bins
k = 3*popN;                                             % # of parmeters estimated
ll = -T/2 * ((log(2*pi) + 1) + log(RSS/T));             % Log-liklihood
AIC = -2*ll + 2*k;                                      % Akaike Information Criterion

% Find first local minimum of AIC score
for ii = 2:length(AIC)
    dAIC(ii-1) = AIC(ii) - AIC(ii-1);
end

% Use local minimum if present, otherwise, use lowest score achieved
if any(dAIC>0)
    locMin = find(dAIC>0, 1);
else
    locMin = find(AIC == min(AIC));
end

% Get parameters from best fit
[subMean, subVar, curveVals] = gaussStats(fitresult{locMin},binCenter);
% subMean = flip(subMean);
% subVar = flip(subVar);
% curveVals = flip(curveVals);

fwhm = 2*sqrt(2*log(2))*sqrt(subVar);

% Make output figures
fig = figure; hold on
fig.Position = [1921 221 1600 783];

if exist("pltTitle", 'var')
    title(pltTitle, 'FontSize', 18, 'Interpreter', 'latex')
else
    title('Data Histogram and Best Gaussian Fit', 'FontSize', 18, 'Interpreter', 'latex')
end

if exist("xaxTitle", 'var')
    xlabel(xaxTitle, 'FontSize', 14, 'Interpreter', 'latex')
else
    xlabel('Bin', 'FontSize', 14, 'Interpreter', 'latex')
end

ylabel('Counts', 'FontSize', 14, 'Interpreter', 'latex')
sumCurve = sum(curveVals,1);
plot(binCenter, sumCurve)
b = bar(binCenter, histOut, [0.5, 0.5, 0.5], 'hist');
b.FaceColor = 'k';
b.FaceAlpha = 0.5;
ax = gca;
ax.Position = [0.05 0.11 0.9 0.815];
ax.XLim = [0 ax.XLim(2)];
ax.YLim = [0 1];
texN = text(0.01+ax.XLim(2),0.85*ax.YLim(2), '$\underline{Fit Results}$', 'FontSize', 16, 'Interpreter', 'latex');
if any(dAIC>0)
    texLab = text(0.01+ax.XLim(2),0.81*ax.YLim(2), ['Gaussians = ', num2str(locMin), '   $R^2 = $', num2str(rsquare(locMin))], 'FontSize', 14, 'Interpreter', 'latex');
else
    texLab = text(0.01+ax.XLim(2),0.81*ax.YLim(2), ['Gaussians = ', num2str(locMin), '*   $R^2 = $', num2str(rsquare(locMin))], 'FontSize', 14, 'Interpreter', 'latex');
end

% Find statistically significantly different peaks (using
% simplified z-statistic and two-tailed p)
for ii = 1:length(subMean)
    for jj = 1:length(subMean)
        z(ii,jj) = (subMean(ii) - subMean(jj)) / sqrt(subVar(ii) + subVar(jj));
        p(ii,jj) = 2*normcdf(z(ii,jj));
    end
end
[peakA, peakB] = find(p<0.05);

% Plot significantly different combination(s)
colorWheel = {'r', 'g', 'b', 'c', 'm', 'y'};
lyneStile = {':', '-.', '--', '-'};
lw = .75;
for ii = 1:length(peakA)
    lw = lw+0.25;
    plt = plot(binCenter, curveVals(peakA(ii),:), [colorWheel{mod(ii-1,6)+1}, lyneStile{mod(ii-1,4)+1}], binCenter, curveVals(peakB(ii),:), [colorWheel{mod(ii-1,6)+1}, lyneStile{mod(ii-1,4)+1}]);
    plt(1).LineWidth = lw;
    plt(2).LineWidth = lw;

    % Plot the means of sigdiff combo(s)
    xLn = xline([subMean(peakA(ii)), subMean(peakB(ii))], ['k', lyneStile{mod(ii-1,4)+1}]);
    xLn(1).LineWidth = 1;
    xLn(2).LineWidth = 1;

    % Find point to center FWHM line around
    dist = abs(binCenter - subMean(peakA(ii)));
    idx = find(dist == min(dist));
    hm(1) = curveVals(peakA(ii),idx)/2;

    dist = abs(binCenter - subMean(peakB(ii)));
    idx = find(dist == min(dist));
    hm(2) = curveVals(peakB(ii),idx)/2;

    % Plot FWHM Line centered over submean
    fw = [subMean(peakA(ii)) - 0.5*fwhm(peakA(ii)), subMean(peakA(ii)) + 0.5*fwhm(peakA(ii))];
    plot(fw, [hm(1), hm(1)], ':k', 'LineWidth', 1)

    fw = [subMean(peakB(ii)) - 0.5*fwhm(peakB(ii)), subMean(peakB(ii)) + 0.5*fwhm(peakB(ii))];
    plot(fw, [hm(2), hm(2)], ':k', 'LineWidth', 1)
end

% Add labels for all sub peaks fit
for ii = 1:locMin
    yInc(ii) = (ii-1)*0.1;
    
    muString =  ['$\mu_', num2str(ii-1) ' = ', num2str(subMean(ii)), '$'];
    texMu = text(0.01+ax.XLim(2),(0.755-yInc(ii))*ax.YLim(2), muString, 'FontSize', 14, 'Interpreter', 'latex');
    
    siString =  ['$\sigma^{2}_{' num2str(ii-1), '} = ', num2str(subVar(ii)), '$'];
    texSi = text(0.01+ax.XLim(2),(0.72-yInc(ii))*ax.YLim(2), siString, 'FontSize', 14, 'Interpreter', 'latex');
end

% Mark significantly different pairs
mark = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j'};
for ii = 1:length(peakA)
    xInc = ii*(.008*ax.XLim(2));
    texSg = text(ax.XLim(2)*1.15 + xInc,(0.7385-yInc(peakA(ii)))*ax.YLim(2), mark{mod(ii-1,10)+1}, 'FontSize', 10, 'Interpreter', 'latex');
end

for ii = 1:length(peakB)
    xInc = ii*(.008*ax.XLim(2));
    texSg = text(ax.XLim(2)*1.15 + xInc,(0.7385-yInc(peakB(ii)))*ax.YLim(2), mark{mod(ii-1,10)+1}, 'FontSize', 10, 'Interpreter', 'latex');
end

% Annotate if AIC had no local minimum
if ~any(dAIC>0)
    texNot = text(0, -0.1*ax.YLim(2), '*No AIC Local Minimum found. Showing curve value for lowest AIC score calculated.','FontSize', 8, 'Interpreter', 'latex');
end

leg = legend('Gaussian Fit', 'Raw Histogram Values', 'FontSize', 16, 'Location', 'northeastoutside', 'Interpreter', 'latex');

hold off;

end