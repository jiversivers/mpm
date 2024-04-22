function flimData = loadFlimFitResults(fname)
% LOADFLIMFITRESULTS loads the FLIM fit results that are exported from
% SPCImage.
%   flimData = LOADFLIMFITRESULTS(fname) returns a struct variable that
%   contains all possible exported results.
%
% Original Code: CAlonzo, loadFlimFitResults.m
% Created: 2014-05-22
% Modified by: Olivia Kolenc, 2020-09-21 -Modified to include wider range of data exported from SPCImage.
% Modified by: Alan Woessner, 11-4-2021 -Created mega function to condense
% multiple functions
% Modified by: Jesse Ivers, 10-10-2022 -Added TIFF support

% Load data from files and place in structure
fileList = dir([fname filesep '*.asc']);
if isempty(fileList)
    fileList = dir([fname, '*.tif*']);
    dtype = '.tif';
else
    dtype = '.asc';
end

for i = 1:size(fileList,1)
    tmpName = strsplit(fileList(i).name,{'_',dtype});
    
    tmpName = tmpName{end-1};
    
    if contains(tmpName,'[%]')
        tmpName = [tmpName(1:end-3),'_per'];
    end
    
    if strcmp(dtype,'.asc')
        flimData.(tmpName) = load([fileList(i).folder,'\',fileList(i).name], '-ascii');
    else
        flimData.(tmpName) = double(imread([fileList(i).folder,'\',fileList(i).name]));
    end
end
% Leave normalization/tauM to examples or demo
% Normalization: Normalize 'a' values based on sum of all 'a' values
%   If only a1, then normalize by max/min values
% TauM: sum(a# (normalized) * t#)
end
