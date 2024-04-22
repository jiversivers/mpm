function flimNameCleaner(directory, ext)

if ~exist('ext', 'var')
    ext = '.asc';
end

% This allows for the ext to be entered with or without the dot and still
% work
x = strsplit(ext, '.');
ext = x{end};

% Get names
fileNames = {dir([directory filesep '*.' ext]).name};

% Create array of clean neams
cleanNames = cellfun(@(x) x{end}, ...
                cellfun(@(x) strsplit(x, '_'), fileNames, 'UniformOutput',false), ...
             'UniformOutput', false);

for ii = 1:numel(fileNames)
    movefile([directory filesep fileNames{ii}], [directory filesep cleanNames{ii}]);
end