function fitData = loadFitData(flimDir)

%{
This data takes the image directory path (containgin .asc or.tiff outputs
from SPCImage for after resorting) as an argument and will output a
structure contiaing all of the fit parameters found. These parameters can
be accessed by dot-indexing with the parameter name.
%}

files = dir(flimDir);
files = files(~[files.isdir]);

[~, param, dtype] = fileparts({files.name});
param = strrep(param, '[%]', 'percent');
assert(all(strcmp(dtype, dtype{1})), 'ERROR LOADING: Ambiguous data set. Data are different types.')

switch dtype{1}
    case '.asc'
        for ii = 1:numel(files)
            try
                fitData.(param{ii}) = load([files(ii).folder filesep files(ii).name], '-ascii');
            end
        end
    case {'tif', 'tiff'}
        for ii = 1:numel(files)
            try
                fitData(params{ii}) = tiffreadVolume([files(ii).folder filesep files(ii).name]);
            end
        end
end
fitData.a1 = abs(fitData.a1);
fitData.a1percent = abs(fitData.a1percent);
fitData.a2 = abs(fitData.a2);
fitData.a2percent = abs(fitData.a2percent);
fitData = flimCombination(fitData);