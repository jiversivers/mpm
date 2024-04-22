clearvars; close all;
errors = cell(0);
%% 
% Simulation parameters
n = 10000;      % Number of iterations for each population count
popN = 1:5;     % Population counts to simulate
fitFlags = {'-r2', '-aic', '-daic', '-ddaic'};  % Selection criteria
fitOpt = listdlg('PromptString', {'What fit criterion would you like to test?'}, 'SelectionMode', 'multiple',...
    'ListString', fitFlags, 'ListSize', [350 100]);


% Prepare for plot
cW = {'r', 'g', 'b', 'c', 'm'};
fig = figure;
fig.Position = [1921 221 1600 783];

%% 
t1 = tic;
% Test each fit criteria
for kk = 1:length(fitOpt)
    % Iterate through each population simulation
    for jj = popN

        % Initiate output data arrays
        AverageofData = zeros([n, 1]);
        StDofData = zeros([n, 1]);
        ActualCount = zeros([n, 1]);
        PredictedCount = zeros([n, 1]);

        % Simulate n data sets 
        for ii = 1:n
            try
                %% Validate HeteroAnaylisis by creating random histograms and fitting. 
                % Create random distribution
                samples = sampleRandomGMM(jj);
                
                % Scale vals to match ORR bounds
                samples = samples/max(samples);
                
                % Calculate histogram
                [N, EDGES] = histcounts(samples, 50);
                N = N(2:end-1);
                centers = 0.5*(EDGES(2:end)+EDGES(1:end-1));
                centers = centers(2:end-1);
                
                % Input into curve fitter
                % Fit GMM
                [curveVals, fitStats, R2] = freqCurve(centers, N, '-i', fitFlags{fitOpt(kk)});
                Count = size(fitStats, 1);
                
                % Append results
                AverageofData(ii*jj) = mean(samples);
                StDofData(ii*jj) = std(samples);
                ActualCount(ii*jj) = jj;
                PredictedCount(ii*jj) = Count;
        
                t = toc(t1);
                estimate = dispTimeRemaining(t, ii*jj*kk, n*length(popN)*length(fitOpt));
                clc
                disp(estimate);
            catch ME
                try
                    warning(getReport(ME))
                    body = ['Error during epoch ' num2str(jj), ' populations testing ', fitFlags{fitOpt(kk)}, ' fitting. ', estimate, 'Error report: ', getReport(ME)];
                    sendolmail('jdivers@uark.edu', 'RE: Simulation Update', body)
                    errors{end+1} = ME;
                catch
                    errors{end+1} = ME;
                end
            end
        end
        
        %% Save Results
        simResults = table(AverageofData, StDofData, ActualCount, PredictedCount);
        writetable(simResults, ['simResults', fitFlags{fitOpt(kk)}, '.xlsx'])


        % Plot each epoch(n)
        scatter3(StDofData, AverageofData, PredictedCount, cW{jj})
        title(['Fits for ', num2str(jj), ' simulated components using ', fitFlags{fitOpt(kk)}])
        ylabel('Average of Simulated Data')
        xlabel('Standard Deviation of Simulated Data')
        zlabel('Predicted Count of Subpopulations')
        ax = gca;
        ax.Position = [0.1300 0.1434 0.7322 0.7816];
        ax.CameraPosition = [-7, -3, 13];
        ax.XLim = [0 1];
        ax.YLim = [0 1];
        ax.ZLim = [1 5];
    
        % Save figure
        try
            exportgraphics(fig, [num2str(jj), fitFlags{fitOpt(kk)}, 'simFig.png'])
        catch
            saveas(fig, [num2str(jj), fitFlags{fitOpt(kk)}, 'simFig.png'])
        end
        
        % Try to send update mail
        try
            body = ['Epoch complete for ' num2str(jj), ' populations testing ', fitFlags{fitOpt(kk)}, ' fitting. ', estimate];
            sendolmail('jdivers@uark.edu', 'RE: Simulation Update', body, {['\\deckard\bmeg\Rajaram-Lab\Ivers,Jesse\UMSCC Radiation Study\Simulation Fit Test\', num2str(jj), fitFlags{kk}, 'simFig.png']})
        catch MEmail
            warning(getReport(MEmail))
        end
    end

end
