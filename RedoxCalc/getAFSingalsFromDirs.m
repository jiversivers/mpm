function varargout = getAFSingalsFromDirs(varargin)

% This funciton will take the two directories (for 755 and 855 images)
% input (in any order) and calculate the redox ouputs from them. If no
% power directory is specified, the program will search for a matching
% powerfile in both image directories and the shared root (if one exists).
% If no power file is found, you will be prompted to select one.
% 
% The first output will be a 1x3 struct with image arrays for each relevent
% signal: NADH, FAD, and SHG. The individual singals can be acessed through
% dot-indexing using those names.
% 
% The second output is the metadata structrue for both images and for the
% power calibraiton file. This consists of the following fields: Laser
% Wavelength  Channels    Gain    Pockels   Date    GridData    ImageSize
% ImageCount  FileNames   Power   These metadata values can be accessed via
% dot indexing into the output structure using these names. The power
% metadata is stored as a map nested beneath Laser and Objective fields.
% The map is accesed using number keys corressponding to the nominal power
% of the measurement.

p = inputParser;
isgoodimg = @(x) all(cellfun(@isfolder, x));
isgoodpwr = @(x) isfolder(x) || isfile(x);
addRequired(p, 'imgDirs', isgoodimg)
addOptional(p, 'pwr', '', isgoodpwr)

p.parse(varargin(1:2), varargin{3:end});

% Get imaging metadata
xmls = cellfun(@(x) [x filesep dir([x filesep '*.xml']).name], p.Results.imgDirs, 'UniformOutput', false);
varargout{2} = cellfun(@readPVxml, xmls);

% Get images
img = cellfun(@imget, p.Results.imgDirs, 'UniformOutput', false);

for ii = [1 2]
    % Get power metadata
    if isempty(p.Results.pwr)
        s = 0;
        searchDirs = horzcat(p.Results.imgDirs{1:2}, {findRoot(p.Results.imgDirs(1:2))});
        while  ~exist('pwrFile', 'var') || isempty(pwrFile)
            s = s+1;
            % Check the s-th searchDir for a pwrFile
            pwrFile = dir([searchDirs{s} filesep '*laserpower*.xlsx']);
            if s == numel(searchDirs)
                [pwrFile, pwrDir] = uigetfile('*.xlsx', 'Select power file.');
                pwrFile = dir([pwrDir fielsep pwrFile]);
            end
        end

    elseif isfolder(p.Results.pwr)
            varargout{2}(ii).Power = mapPower(getPower(p.Results.pwr, varargout{2}(ii).Date), varargout{2}(ii).Wavelength);
    else
            error('Unrecognized power input')
    end

    % Normalize img
    if ismissing(varargout{2}(ii).Power) || isempty(varargout{2}(ii).Power)
        warning('Could not normalize image without power data.')
        img{ii} = NaN(size(img{ii}));
    else
        img{ii} = ImageNormalize(img{ii}, varargout{2}(ii), 'fluorescein');
    end
end

% Figure out which is which
varargout{1}.NADH = img{[varargout{2}.Wavelength] == 755}(:,:,3); %NADH
varargout{1}.FAD = img{[varargout{2}.Wavelength] == 855}(:,:,2); %FAD
varargout{1}.SHG = img{[varargout{2}.Wavelength] == 855}(:,:,4); %SHG
varargout{1}.all = cat(3, img{[varargout{2}.Wavelength] == 755}, img{[varargout{2}.Wavelength] == 855});
