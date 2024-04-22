function pwr = getPower(power, dateTaken)

if isfolder(power)

% If no date input, assume power file is only one in directory and is named
% something with power in the name
    if ~exist('dateTaken', 'var')
        name = dir([power, filesep, '*power*.xlsx']).name;
        pwrpath = dir([power, filesep, '*power*.xlsx']).folder;

% Else, find the power file by its date 
    else 
        ymdDate = num2str(yyyymmdd(dateTaken(1)));
        
        % Reformat date to amtch save name format
        if strcmp(ymdDate(5), '0')
            mdyDate = [ymdDate(6), ymdDate(7:8), ymdDate(1:4)]; % Disregard leading zero
        else
            mdyDate = [ymdDate(5:6), ymdDate(7:8), ymdDate(1:4)];
        end
    
        try
            name = dir([power, filesep, '*',mdyDate,'*.xlsx']).name; 
            pwrpath = dir([power, filesep, '*',mdyDate,'*.xlsx']).folder;
        catch
            warning('Power file for date not found.')
            pwr = missing;
            return
        end
    end

    % Read power data from 
    sheets = sheetnames([pwrpath, filesep, name]);
    if numel(sheets) > 1
        pwr = readcell([pwrpath, filesep, name], 'Sheet', sheets{2});
    else
        pwr = readcell([pwrpath, filesep, name]);
    end

elseif isfile(power)
    sheets = sheetnames(power);
    if numel(sheets) > 1
        pwr=readcell(power, 'Sheet', sheets{2});
    else
        pwr = readcell(power);
    end

else
    error('Power input must be a file or a directory of containing a power file.')
end