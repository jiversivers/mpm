clearvars; close all;

%% HeteroAnalysis
% This program will load and analyze in batches previously processed ORR
% data (from batchORR and/or related functions). Each sample is
% automatically extracted and run through weightedHistFits.
% weightedHistFits creates an intensity map (average of NADH and FAD
% fluroescence) weighted histogram of ORR values in an image. This
% histogram is then fit with a 1-D gaussian mixture model (superposition of
% Gaussian curves). The number of components in the fit is determined by
% the first positive inflection point of the AIC score (i.e. the first
% point where the rate of AIC-decrease slows). The raw histogram and the
% fit are both plotted automattically, along with the individual Gaussians.
% All details regarding the fit will be included in the figure (mean,
% variance, and portion of each curve). along with an annotation of the
% number of components in the fit. All of this data (the inidividual curve
% stats) are saved for each FOV. Further, a sturct will be created with  and the figure are saved to the
% sample folder.

% Updated: 28-Sep-2022 by Jesse Ivers

% Updated with FLIM options


%% Get analysis procedure inputs

imgOpt = questdlg('What image type are you analyzing?', ...
    'Image Type', ...
    'Redox maps', 'Lifetime maps', 'Redox maps');

choice = {'R^2 > 0.95', 'Minimum AICc', 'Maximum dAICc (greatest decrease in improvment)', ...
    'First ddAICc > 0 (first inflection point)'}

fitOpt = listdlg('PromptString', {'What fit criterion would you like to use?'}, 'SelectionMode', 'single',...
    'ListString', choice, 'ListSize', [350 100]);

if fitOpt == 1
    fitFlag = '-r2'
elseif fitOpt == 2
    fitFlag = '-aic';
elseif fitOpt == 3
    fitFlag = '-daic';
else
    fitFlag = '-ddaic';
end

% % Get starting directory to process out of
disp('Navigate to processed root data directory')
procDir = uigetdir(pwd,'Select root directory of processed data');
cd(procDir)
clc

disp('Locating directories of data to be analyzed...')

if strcmp(imgOpt,'Redox maps')
    mapName = '*_ORRmap.mat';
    intName = '_data_struct.mat';
else
    a1mapName = '*_a1[%].tiff';
    a2mapName = '*_a2[%].tiff';
    t1mapName = '*_t1.tiff';
    t2mapName = '*_t2.tiff';
    intName = '_intensity_image.tif';
end

proc_terms = FolderFinder(intName);

disp('Creating figures and saving outputs...')

%%

% Prepare compiled output table variables
CellLine = cell(length(proc_terms),1);
Group = cell(length(proc_terms),1);
AnimalID_Sample = cell(length(proc_terms),1);
FOV = cell(length(proc_terms),1);
Count = cell(length(proc_terms),1);
GOF = cell(length(proc_terms),1);
R2 = cell(length(proc_terms),1);
curveStats = cell(length(proc_terms),1);

t = 0;
tStart = tic;
errors = cell(0);
for jj = 1:length(proc_terms)
    cd(proc_terms{jj})    
    try 
        if strcmp(imgOpt, 'Redox maps')
    
            % Load variables
            load(dir(mapName).name);
            load(dir(['*', intName]).name);
            
            % Generate intensity map (avg of NADH and FAD intensities)
            for ii = 1:length(Data)
                if Data(ii).Wavelength == '755'
                    nadh = Data(ii).UseImg(:,:,3);
                elseif Data(ii).Wavelength == '855'
                    fad = Data(ii).UseImg(:,:,2);
                end
            end
            
            int = (nadh+fad)./2;
        
            % Prepare histogram data
            [N, centers] = weightedHistogram(map, int, 50);
            N = N(2:end-1);
            N = N-min(N);
            N = N/sum(N);
            centers = centers(2:end-1);
    
            % Fit GMM
            [curveVals, fitStats, R2{jj}] = freqCurve(centers, N, '-i', fitFlag);
            Count{jj} = size(fitStats, 1);
    
            % Make fiugre
            fig = freqCurvePlot([centers; N]', curveVals, fitStats, R2{jj}, 'Intensity-weighted Histogram of ORR', 'ORR Value', 'Weighted Frequency');
    
            % Creating file names
            name = pathpieces(proc_terms{jj});
            name = matlab.lang.makeValidName(name);
            
            % Save inidividual curve stats
            curveStats{jj} = struct('subMean', fitStats(:,1), 'subVariance', fitStats(:,2), 'Portion', fitStats(:,3));
            if exist([name{end}, '_curveStats.xlsx'],"file")
                delete([name{end}, '_curveStats.xlsx'])
            end
            writetable(struct2table(curveStats{jj}), [name{end}, '_curveStats.xlsx'])
            
            % Save figure
            exportgraphics(fig, [name{end}, '_HistFig','.png'])
            
            % Append count labels to compiled data table
            CellLine{jj} = name{end-3};
            Group{jj} = name{end-2};
            AnimalID_Sample{jj} = name{end-1};
            FOV{jj} = name{end};
            
            clc
            dispTimeRemaining(t,jj,length(proc_terms))
            disp(['Succesfully saved figure for ', name{end-1}, ' ', name{end}, '.'])
            close all; clear fig wHistCounts curveVals fitStats

        
        else
            % Load variables
            int = double(imread(dir(['*', intName, '*']).name));
            a1 = double(imread(dir(a1mapName).name));
            a2 = double(imread(dir(a2mapName).name));
            t1 = double(imread(dir(t1mapName).name));
            t2 = double(imread(dir(t2mapName).name));

            % Create Tau_m and a1/a2 maps
            map{1} = a1;
            map{2} = a1.*t1 + a2.*t2;

            % Create labels
            label{1,1} = 'Intensity Weighted Histogram of $$ \alpha_1 $$';
            label{1,2} = 'Intensity Weighted Histogram of $$ \tau_m $$';
            label{2,1} = '$$ \alpha_1 $$';
            label{2,2} = '$$ \tau_m (ps) $$';
            label{3,1} = 'Total Weight (a.u.)';
            label{3,2} = 'Total Weight (a.u.)';
            label{4,1} = 'a1';
            label{4,2} = 'taum';

            for kk = 1:2
                % Prepare histogram data
                [N, centers] = weightedHistogram(map{kk}, int, 50);
                N = N(2:end-1);
                N = N-min(N);
                N = N/sum(N);
                centers = centers(2:end-1);
        
                % Fit GMM
                [curveVals, fitStats, R2{jj,kk}] = freqCurve(centers, N, '-f', fitFlag);
                Count{jj,kk} = size(fitStats, 1);
        
                % Make fiugre
                fig = freqCurvePlot([centers; N]', curveVals, fitStats, R2{jj,kk}, label{1,kk}, label{2,kk}, label{3,kk});
        
                % Creating file names
                name = pathpieces(proc_terms{jj});
                name = matlab.lang.makeValidName(name);
                
                % Save inidividual curve stats
                curveStats{jj,kk} = struct('subMean', fitStats(:,1), 'subVariance', fitStats(:,2), 'Portion', fitStats(:,3));
                if exist([name{end}, '_curveStats.xlsx'],"file")
                    delete([name{end}, '_curveStats.xlsx'])
                end
                writetable(struct2table(curveStats{jj,kk}), [name{end}, label{4,kk}, '_curveStats.xlsx'])
                
                % Save figure
                exportgraphics(fig, [name{end}, label{4,kk}, '_HistFig','.png'])
                
                % Append count labels to compiled data table
                CellLine{jj} = name{end-3};
                Group{jj} = name{end-2};
                AnimalID_Sample{jj} = name{end-1};
                FOV{jj} = name{end};
                
                clc
                dispTimeRemaining(t,jj,length(proc_terms))
                disp(['Succesfully saved ', label{4,kk}, ' figure for ', name{end}, '.'])
                close all; clear fig wHistCounts curveVals fitStats
            end
        end

    catch ME
        clc
        disp('Failed for some reason...see below:')
        warning(getReport(ME));
        errors{end+1} = ME;
    end

    t = toc(tStart);

end

% Save compiled  table
outTab = table(CellLine, Group, AnimalID_Sample, FOV, Count, R2);
writetable(outTab,[procDir, filesep, 'SPA_Counts.xlsx']);

disp('Output table saved.')
disp('Done')



    