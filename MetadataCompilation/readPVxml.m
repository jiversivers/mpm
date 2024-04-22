function pvData = readPVxml(filename)

if ~isfile(filename) && isfolder(filename)
    filename = [filename filesep dir([filename filesep '*.xml']).name];
end

% Parse xml into nested struct
xml = xml2struct_nonnative(filename);

% Extract key names from all the nested metadata structs
keys = cellfun(@(x) x.Attributes.key, xml.PVScan.PVStateShard.PVStateValue, 'UniformOutput', false);

% Pluck out nested struct to make indexing in cleaner
metaData = xml.PVScan.PVStateShard.PVStateValue;

%% Imaging Modality
mode = metaData{strcmp(keys, 'activeMode')}.Attributes.value;

%% FLIM Specifics
if strcmpi(mode, 'FLIM')
    try
        dt = str2double(xml.PVScan.AcquisitionData.Attributes.FLIMTimeBinNanoseconds);
    catch
        dt = 'missing';
    end
else
    dt = 'N/A';
end


%% Wavelength
lambda = str2double(metaData{strcmp(keys, 'laserWavelength')}.IndexedValue.Attributes.value);

%% Power (Pockels)
pock = str2double(metaData{strcmp(keys, 'laserPower')}.IndexedValue.Attributes.value);

%% PMT Info
% Color
color = cellfun(@(x) x.Attributes.description, metaData{strcmp(keys, 'pmtGain')}.IndexedValue, 'UniformOutput', false);

% Gain
pmt = cellfun(@(x) str2double(x.Attributes.value), metaData{strcmp(keys, 'pmtGain')}.IndexedValue);

%% Laser type
laser = metaData{strcmp(keys, 'laserWavelength')}.IndexedValue.Attributes.description;
if ~strcmp(laser, 'InsightX3')
    laser = 'MaiTai';
end

%% Image Count
imgCount = numel(xml.PVScan.Sequence);

%% Image Size
imgSize = [str2double(metaData{strcmp(keys, 'pixelsPerLine')}.Attributes.value)  ...
           str2double(metaData{strcmp(keys, 'linesPerFrame')}.Attributes.value)];

%% Image scale
imgScale = [str2double(metaData{strcmp(keys, 'micronsPerPixel')}.IndexedValue{1}.Attributes.value) ...
            str2double(metaData{strcmp(keys, 'micronsPerPixel')}.IndexedValue{2}.Attributes.value)];

%% Grid Data
seq = xml.PVScan.Sequence;
if imgCount > 1
    try
        gridData(1) = str2double(seq{1,1}.Attributes.xYStageGridNumXPositions);
        gridData(2) = str2double(seq{1,1}.Attributes.xYStageGridNumYPositions);
        gridData(3) = str2double(seq{1,1}.Attributes.xYStageGridOverlapPercentage);
    catch
        gridData = 'Undefined';
    end
else
    gridData = [1 1 100];
end

%% Filenames
filenames = cell(imgCount, numel(color));
if imgCount > 1
        for cyc = 1:imgCount
            for ch = 1:numel(color)
                try
                    filenames{cyc, ch} = seq{1,cyc}.Frame.File{1,ch}.Attributes.filename;
                catch
                    filenames{cyc, ch} = 'Undefined';
                end
            end
        end
else
    filenames = cellfun(@(x) x.Attributes.filename, seq.Frame.File, 'UniformOutput', false);
end

%% Date of imaging
dateTaken = datetime(xml.PVScan.Attributes.date,'InputFormat','MM/dd/uuuu hh:mm:ss aa'); 

%% Final struct
pvData = struct('Mode', mode, ...
                'Laser', laser, ...
                'Wavelength', lambda, ...
                'Channels', {color}, ...
                'Gain', pmt, ...
                'Pockels', pock, ...
                'Date', dateTaken, ...
                'GridData', gridData, ...
                'ImageSize', imgSize, ...
                'ImageScale', imgScale, ...
                'TemporalBinWidth', dt,...
                'ImageCount', imgCount, ...
                'FileNames', {filenames});

