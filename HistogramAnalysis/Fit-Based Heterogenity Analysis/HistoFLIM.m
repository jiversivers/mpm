clearvars; close all;

%% HistoFLIM
% This program will load and analyze in batches previously processed FLIM
% data (from SPCImage and postSPCImage or flimAscii2Tiff). Each sample is
% automatically extracted and run through weightedHistFits.
% weightedHistFits creates an intensity-map-weighted histogram of lifetime
% values in an image. This histogram is then fit with Gaussian curve(s).
% The raw histogram and the fit are both plotted automattically. The
% minimum number of superimposed Gaussian curves is used such that the R
% squared value of the fit is equal to or greater than 95% (0.95). This fit
% is then deconstructed into its individual Gaussian curves which are then
% compared through a simplified z-test. This z-test takes the square
% difference of the means divided by the square root of the sum of the
% variances. From this z score, a p value is determined. For those pairs of
% individual curves with a p value of no more than 0.05, the associated
% individual curves are added to the plot, along with markers for their
% corresponding means and FWHM values. Alongside this plot, mean and
% variance for all Gaussian curves in the model are listed, along with an
% annotation of the number of components in the fit. All of this data (the
% inidividual curve stats as well as comparison results) and the figure are
% saved to the sample folder.

% Updated: 26-Jul-2022 by Jesse Ivers

% Switched to freqCurve in place of WeigthedHistFits. Added estimated time
% remaining. Added error capturing structure. 


%% 

% % Get starting directory to process out of
disp('Navigate to processed root data directory')
procDir = uigetdir(pwd,'Select top-most processed data directory in directory tree to be analyzed');
clc

disp('Navigate to location for save directory')
saveDir = uigetdir(procDir,'Select save location for output directory');
cd(saveDir)
if ~exist('HistoFLIMOut', 'dir')
    mkdir('HistoFLIMOut')
end
cd(procDir)
clc

disp('Locating directories of data to be analyzed...')
proc_terms = FolderFinder('_intensity_image.tif');

totTime = 0;
%% Send each sample through WeightHistFits function
disp('Calculating fit parameters and creating figures...')
tic;
errors = cell(0);
for jj = 1:length(proc_terms)
    cd(proc_terms{jj})
    
    try 
        % Lifetime variable names (row 1 for plot, row 2 for axis, row 3 for file)
        names = {'$\tau_1$ (ps)', '$\tau_2$ (ps)', '$\alpha_1$ (\%)','$\alpha_2$ (\%)', '$\tau_m$ (ps)', '$\frac{\alpha_1}{\alpha_2}$'; ...
                'Lifetime',         'Lifetime',     'Fraction',     'Fraction',         'Lifetime',     'Relative Fraction';...
                't1',               't2',           'a1[%]',        'a2[%]',            'tm',           'a1a2ratio'};
        % Load variables
        par{1} = double(imread(dir('*_t1.tif*').name));
        par{2} = double(imread(dir('*_t2.tif*').name));
        par{3} = double(imread(dir('*_a1[%].tif*').name));
        par{4} = double(imread(dir('*_a2[%].tif*').name));
        par{7} = double(imread(dir('*_intensity_image.tif*').name));
        
        % Calculate mean and ratio
        par{5} = (par{3}/100).*par{1} + (par{4}/100).*par{2};
        par{6} = par{3}./par{4};
        
        % Each row is fit for different var, each column is different
        % output (figure, submean, subvar, zscore, pvalue, curvestats)
        fitData = cell(6);
        
        % Workhorse function...see code for more details
        for ii = 1:6
            [fitData{ii,1}, fitData{ii,2}, fitData{ii,3}, fitData{ii,4}, fitData{ii,5}] = freqCurve(par{ii}, par{7}, ['Intensity-weighted Histogram and Gaussian fit of ', names{1,ii}], names{2,ii});
            
            % Save all
            fitData{ii,6} = struct('subMean', fitData{ii,2}', 'subVariance', fitData{ii,3}', ...
                                'zScore', fitData{ii,4}, 'pValue', fitData{ii,5});
            writetable(struct2table(fitData{ii,6}), [saveDir, filesep, 'HistoFLIMOut', filesep, proc_terms{jj}(length(procDir)+2:end), '_curveStats_', names{3,ii}, '.xlsx'])
            exportgraphics(fitData{ii,1}, [saveDir, filesep, 'HistoFLIMOut', filesep, proc_terms{jj}(length(procDir)+2:end), '_HistFig_', names{3,ii},'.png'])
        end

        save([saveDir, filesep, 'HistoFLIMOut', filesep, proc_terms{jj}(length(procDir)+2:end),'fitData.mat'],'fitData')
        clc
        disp(['Succesfully saved all for ', proc_terms{jj}(length(procDir)+2:end), '.'])

    catch ME
        clc
        disp('Failed for some reason...see below:')
        warning(getReport(ME));
        errors{end+1} = ME;
    end

    t = toc;
    totTime = totTime + t;
    remTime = round((totTime / jj) * (length(proc_terms) - jj));
    scd = mod(remTime,60); 
    min = mod(((remTime-scd)/60),60);
    hrs = (remTime-scd-(60*min))/3600;
    if hrs ~=0
        disTime = [num2str(hrs),' hours, ', num2str(min), ' minutes remaining.'];
    else
        disTime = [num2str(min), ' minutes, ', num2str(scd), ' seconds remaining.'];
    end
    
    disp([num2str(jj), '/', num2str(length(proc_terms)), ' completed.']);
    disp(['Est. time remaining: ', disTime])
    close all;
    clear fitData par;
end
cd(saveDir)

if exist("errors",'var') && numel(errors) > 0
    disp('Done, but the following errors occured:')
    for ii = 1:numel(errors)
        warning(getReport(errors{ii}))
    end
else
    disp('Done.')
end
