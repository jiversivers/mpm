function fig = freqCurvePlot(histCounts, curveVals, fitStats, R2, pltTitle, xAxLabel, yAxLabel)

%% freqCurvePlot 
% This function takes data from the freqCurve analysis and plots the full
% GMM fit along with the individual component curves. Histogram data will
% also be included as a scatter underlayed with the curves. The histCounts
% data is input as a n x 2 array where the first column is the bin value
% and the second column is the corresponding frequency. Alongside this
% plot, mean variance, and portion for all Gaussian curves in the model are
% listed, along with an annotation of the number of components in the fit.

% Last Updated: 11-Aug-2022 by Jesse Ivers
% Reverted back to histogram style for raw data.

%%

% Parse input data
binVal = histCounts(:,1);
freq = histCounts(:,2);
n = min(size(curveVals));
mu = fitStats(:,1);
sig = fitStats(:,2);
port = fitStats(:,3);

% Init & size figures
fig = figure; hold on
fig.Position = [1921 221 1600 783];

% Label plot
% Varargin Title
if exist("pltTitle", 'var')
    title(pltTitle, 'FontSize', 18, 'Interpreter', 'latex')
else
    title('Data Histogram and Best Gaussian Fit', 'FontSize', 18, 'Interpreter', 'latex')
end

% Vargargin x Axis label
if exist("xAxLabel", 'var')
    xlabel(xAxLabel, 'FontSize', 14, 'Interpreter', 'latex')
else
    xlabel('Bin', 'FontSize', 14, 'Interpreter', 'latex')
end

% Vargargin y Axis label
if exist("yAxLabel", 'var')
    ylabel(yAxLabel, 'FontSize', 14, 'Interpreter', 'latex')
else
    ylabel('Frequency', 'FontSize', 14, 'Interpreter', 'latex')
end

% Plot histogram data
br = bar(binVal, freq,'hist');
br.FaceColor = [0.75 0.75 0.75];
br.EdgeColor = [0.5 0.5 0.5];
br.FaceAlpha = 0.5;
br.EdgeAlpha = 0.5;

% Plot summed fit
sumCurve = sum(curveVals,1);
plot(binVal, sumCurve)

% Add description of fit to figure
ax = gca;
ax.Position = [0.05 0.11 0.9 0.815];
ax.XLim = [0 ax.XLim(2)];
ax.YLim = [0 max(histCounts(:,2))];
text(0.01+ax.XLim(2),.98*ax.YLim(2), '$\underline{Fit Results}$', 'FontSize', 16, 'Interpreter', 'latex');
text(0.01+ax.XLim(2), 0.94*ax.YLim(2), ...
    ['Components = ', num2str(n), ' $\enspace R^2 = ', num2str(R2), '$'], ...
    'FontSize', 14, 'Interpreter', 'latex');

% Plot and list component fits
if n>1
    curveNames = cell(1,n);
    for ii = 1:n
        % Plot
        plot(binVal, curveVals(ii,:));
        curveNames{ii} = ['Component ', num2str(ii)];
    
        % Labels
        yInc = (ii-1)*0.15;
        
        muString =  ['$\mu_', num2str(ii) ' = ', num2str(mu(ii)), '$'];
        text(0.01+ax.XLim(2),(0.885-yInc)*ax.YLim(2), muString, 'FontSize', 14, 'Interpreter', 'latex');
        
        siString =  ['$\sigma^{2}_{' num2str(ii), '} = ', num2str(sig(ii)), '$'];
        text(0.01+ax.XLim(2),(0.85-yInc)*ax.YLim(2), siString, 'FontSize', 14, 'Interpreter', 'latex');
    
        piString = ['$\pi_{' num2str(ii), '} = ', num2str(round(port(ii),2)*100), '\%$'];
        text(0.01+ax.XLim(2),(0.815-yInc)*ax.YLim(2), piString, 'FontSize', 14, 'Interpreter', 'latex');
end
    legend(['Raw Histogram Values', 'Gaussian Fit', curveNames], 'FontSize', 16, 'Location', 'southeastoutside', 'Interpreter', 'latex');
else
    muString =  ['$\mu = ', num2str(mu), '$'];
    text(0.01+ax.XLim(2),0.885*ax.YLim(2), muString, 'FontSize', 14, 'Interpreter', 'latex');
    
    siString =  ['$\sigma^{2} = ', num2str(sig), '$'];
    text(0.01+ax.XLim(2),0.85*ax.YLim(2), siString, 'FontSize', 14, 'Interpreter', 'latex');

    piString = ['$\pi_ = ', num2str(round(port,2)*100), '\%$'];
    text(0.01+ax.XLim(2),0.815*ax.YLim(2), piString, 'FontSize', 14, 'Interpreter', 'latex');
    legend(['Raw Histogram Values', 'Gaussian Fit'], 'FontSize', 16, 'Location', 'southeastoutside', 'Interpreter', 'latex');
end

hold off;
end

