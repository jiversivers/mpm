function powerMap = mapPower(powerXL, lambda)

if ~iscell(powerXL)
    powerMap = missing;
    return
end

% Converting power data to Matlab structure
n = cellfun(@isnumeric, powerXL);
% Determining row indecies in spreadsheet
op_st_i = find(strcmp(powerXL,'Objective Power')) + find(n(find(strcmp(powerXL,'Objective Power')):end, 2), 1) - 1; % Obj power start row
op_fi_i =size(powerXL, 1); % Obj power end row should be last row of excel
lp_st_i = find(strcmp(powerXL,'Laser Power')) + find(n(find(strcmp(powerXL,'Laser Power')):end, 2), 1) - 1; % Laser power start row
lp_fi_i = lp_st_i + (op_fi_i - op_st_i); % Laser power end row (should have same number of rows as Obj Power)


% Determine column for given lambda
kk = find(sum(cell2mat(cellfun(@(x)isnumeric(x) && x==lambda, powerXL, 'UniformOutput', false))), 1);  

%% Laser Power
key = powerXL(lp_st_i:lp_fi_i,kk);
pwr = powerXL(lp_st_i:lp_fi_i,kk+1);

%Clear missing values
missingPwr = cellfun(@ismissing, powerXL(lp_st_i:lp_fi_i,kk+1));
missingKey = cellfun(@ismissing, powerXL(lp_st_i:lp_fi_i,kk));
key(missingPwr|missingKey) = [];
pwr(missingPwr|missingKey) = [];

% Map
LP = containers.Map(key, pwr);

%% Obj Power
key = powerXL(op_st_i:op_fi_i,kk);
pwr = powerXL(op_st_i:op_fi_i,kk+1);

%Clear missing values
missingPwr = cellfun(@ismissing, powerXL(lp_st_i:lp_fi_i,kk+1));
missingKey = cellfun(@ismissing, powerXL(lp_st_i:lp_fi_i,kk));
key(missingPwr|missingKey) = [];
pwr(missingPwr|missingKey) = [];

% Map
OP = containers.Map(key, pwr);

% Compile map
powerMap = struct('Laser', LP, 'Objective', OP);
end