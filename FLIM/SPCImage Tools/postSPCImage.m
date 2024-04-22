%% Post-SPCImage FLIM dealer 
% This program wil handle the batch directory that contains the sperate
% sample files after SPCImage exports. You will be asked to select the
% batch directory and where you would like to move the sorted root.You will
% also be asked if you would like to convert ascii files to tiffs. If you
% select yes, the ascii files will be deleted after conversion. All files
% will be sorted back into individual sample folders and renamed with to
% only contain the relavent parameter (as the sample information will be
% preserved in the parent directory of each file).

% Get user input for processing
batchDir = uigetdir(pwd, 'Select batch directory...');
moveDir = uigetdir(batchDir, 'Select sort root...');
copyDecayOpt = questdlg('Do you want to copy decay and fit files to the sorted directory?', ...
     'Copy decay?', ...
     'Yes', 'No', 'Yes');
convertOpt = questdlg('Do you want to convert ASCII files to TIFF?', ...
    'Convert?', ...
    'Yes', 'No', 'No');

% Convert all ASCII to TIFF
if strcmp(convertOpt, 'Yes') 
    ascii2tiff(terminals{ii})

    % Check for succesful conversion and delete old ASCII
    % Get .asc file check list
    ascFiles = {dir('*.asc').name};
    [~, ascFileName, ~] = fileparts(ascFiles);
    
    % Get .tiff file cross check list
    tifFiles = {dir('*.tif*').name};
    [~, tifFileName, ~] = fileparts(tifFiles);

    % Make sure each ascii file has a tiff counterpart then del
    for jj = 1:length(ascFileName)
        if sum(strcmp(ascFileName{jj}, tifFileName)) == 1
            delete(ascFiles{jj})
        elseif sum(strcmp(ascFileName{jj}, tifFileName)) > 1
            warning(['Multiple .tiff matches found for ', ascFileName{jj}, '. This ascii file will not be deleted.'])
        else  
            warning(['Unable to confirm .tiff conversion file for ', ascFileName{jj}, '. This ascii file will not be deleted.'])
        end
    end
end

% Sort back into seperate directories
% Find the base names in the directory from the oringal .sdt files.
namesSDT = {dir([batchDir filesep '*.sdt']).name};
namesSDT = cellfun(@(x) [batchDir filesep x], namesSDT, 'UniformOutput', false);

% Process through all data sets present
for ii = 1:length(namesSDT)
        
    % Find all the files of the current set
    [~, name, ~] = fileparts(namesSDT{ii});
    searchName = ['*', name, '*'];
    groupNames = {dir([batchDir filesep searchName]).name};
    groupNames = cellfun(@(x) [batchDir filesep x], groupNames, 'UniformOutput', false);
    nameParts =  cellfun(@(x) strsplit(x, '_'), groupNames, 'UniformOutput', false);

    % Move each file to newly created folder and rename w/ the parameter
    for jj = 1:length(groupNames)
        if contains(nameParts{jj}{end}, '.img')
            if strcmp(copyDecayOpt, 'No')
                continue
            end
            paramName = 'fit.img';
        elseif contains(nameParts{jj}{end}, '.sdt')
           if strcmp(copyDecayOpt, 'No')
                continue
            end
            paramName = 'decay.sdt';
        else
            paramName = nameParts{jj}{end};
        end

        % Create new directory for current set
        destination = [moveDir filesep namemaker(nameParts{jj}(2:end), filesep)];
        if ~exist(destination, 'dir')
            mkdir(destination);
        end

        movefile(groupNames{jj}, [destination filesep paramName]);
    end
end



