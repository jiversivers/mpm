source = uigetdir('Select root directory of data to be analyzed', pwd);
np = strsplit(source, filesep);
destination = strrep(source, np{end}, 'Processed');
%%%%%%%%%%%%%%%
%% Stitching %%
%%%%%%%%%%%%%%%
%% Open ImageJ through MATLAB
if ~exist('IJM', 'var') || ~isa(IJM, 'net.imagej.matlab.ImageJMATLABCommands')
    ImageJ
end

% Find all the TSeries (regardless of mode)
TSeriesIHC = FolderFinder('References', source);
TSeriesAF = FolderFinder('References', '\\deckard\BMEG\Rajaram-Lab\Ivers,Jesse\UMSCC Radiation Study\MPM\Raw');
TSeries = horzcat(TSeriesIHC, TSeriesAF);
TSeries = TSeries(contains(TSeries, 'TSeries'));
pp = pathpieces(TSeries);

% Back each one out to the sample dir
for ii = 1:size(pp, 1)
    TSeries{ii} = namemaker(pp(ii, 1:sum(~cellfun(@isempty, pp(ii,:)))-1), filesep);
    if strcmp(pp(ii, 1), 'deckard')
        TSeries{ii} = strcat([filesep filesep], TSeries{ii});
    end
end
TSeries = unique(TSeries);

%% Stitch all atlases
errors = {};
for ii = 1:numel(TSeries)
    fprintf('%d', ii);
    ts = TSeries{ii};
    % Get inidviudal atlas folders
    tsDir = FolderFinder('References', ts);
    tsDir =tsDir([contains(tsDir, 'TSeries')]);
    saveDir = regexprep(tsDir, expr, 'D:\\UMSCC Rad Study\\Processed\\$1\\$2\\$3\\$4\\$5\\MPM\\TSeries\\$6');

    % Stitch atlas from metadata
    for jj = 1:numel(tsDir)
        td = tsDir{jj};
        sd = saveDir{jj};
        if ~exist(sd, 'dir')
            mkdir(sd)
        end
        try
            % Stitch the image
            ijXMLStitching(td, sd);
        catch ME
            errors{end+1, 1} = ts;
            errors{end, 2} = ME;
            warning(getReport(ME))
        end
    end
    disp([num2str(jj) ' of ' num2str(numel(tsDir)) ' atlases from sample complete.'])
end
disp([num2str(ii) ' of ' num2str(numel(TSeries)) ' samples complete.'])

%%%%%%%%%%%%%%%%%%%
%% Normalization %%
%%%%%%%%%%%%%%%%%%%
% Get processed folders
Stit = FolderFinder('img_t1_z1_c1', [destination filesep '..' filesep 'Stitched']);
np = cellfun(@(x) strsplit(x, filesep), Stit, 'UniformOutput', false);
np = vertcat(np{cellfun(@(x) size(x, 2) == 11, np)});
%%
for ii = 1:numel(Stit)
%     if exist([Stit{ii} filesep 'normalizedSingle.tiff'], 'file')
%         disp([Stit{ii} ' already normed...skipping'])
%         continue
%     end
    fprintf('Normalizing %d/%d\n', ii, numel(Stit));
    img = loadStitchStack(Stit{ii});
    md = readPVxml([Stit{ii} filesep dir([Stit{ii} filesep '*.xml']).name]);
    md.Power = mapPower(getPower('D:\UMSCC Rad Study\Single Day Power data', md.Date), md.Wavelength);
    img = ImageNormalize(img, md, 'fluorescein');
    TiffSingleWrite(img, [Stit{ii} filesep 'normalizedSingle.tiff'], 'overwrite')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Intramodal Stacking %%
%%%%%%%%%%%%%%%%%%%%%%%%%
norm = FolderFinder('normalizedSingle.tiff', [destination filesep '..' filesep 'Stitched']);
np = pathpieces(norm);
np = vertcat(np{cellfun(@(x) size(x, 2)==11, np)});

% Get image pairs for registration (ID\Mode\Sample)
pairs = unique(cellfun(@(u,v,w,x,y,z) [u filesep v filesep w filesep x filesep y filesep z], ...
    np(:,4), np(:,5), np(:,6), np(:,7), np(:,8), np(:, 9), ...
    'UniformOutput', false));

warning('off', 'MATLAB:imagesci:Tiff:missingExtraSamples')
errors = cell(0, 2);
imgs = cell(2, 1);
for ii = 1:numel(pairs)
%     if exist([destination filesep '..' filesep 'Stitched' filesep pairs{ii} filesep 'TSeries\norm_reg_Stack.tiff'], 'file')
%         continue
%     end
    
    try
        fprintf('Processing %s. Pair %d of %d.\n', pairs{ii}, ii, numel(pairs))
        % Pair images from same sample and modality
        imgpair = cellfun(@(x,y) [x filesep y], ...
            {dir([destination filesep '..' filesep 'Stitched' filesep pairs{ii} filesep 'TSeries' filesep '*5*']).folder}, ...
            {dir([destination filesep '..' filesep 'Stitched' filesep pairs{ii} filesep 'TSeries' filesep '*5*']).name}, ...
            'UniformOutput', false);
        
        % Load pair (755 will always be idx 1, making it the fixed img)
        imgs{1} = tiffreadVolume([imgpair{ [contains(imgpair, '755')]} filesep 'normalizedSingle.tiff']);
        imgs{2} = tiffreadVolume([imgpair{~[contains(imgpair, '755')]} filesep 'normalizedSingle.tiff']);
        
        % Register the pair (only translation for intramodal registration)
        pairOut = threshRegister(imgs{:}, '-tra');
        
        % Save registered stack
        if ~exist([destination filesep pairs{ii} filesep 'TSeries'], 'dir')
            mkdir([destination filesep pairs{ii} filesep 'TSeries'])
        end
        TiffSingleWrite(pairOut, [destination filesep '..' filesep 'Stitched' filesep pairs{ii} filesep 'TSeries\norm_reg_Stack.tiff'], 1)
        disp(['complete' pairs{ii}])
    catch ME
        errors{end+1, 1} = pairs{ii};
        errors{end, 2} = ME;
        disp(['error on ' pairs{ii}])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Manual Masking Tracing %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get registered stack locations
RegStacks = FolderFinder('norm_reg_Stack.tiff', destination);
maskStack = FolderFinder('maskStack.tiff', destination);

RegStacks = RegStacks(~[contains(RegStacks, maskStack)]);

%%
% minibatching to improve time between sets of images. Shuffling to limit
% any systematic error in tracing from working at different times
minibatchsize = 32;
batchIdx = indexBatch(RegStacks, minibatchsize, 'shuffle', true);
impath = cellfun(@(x) RegStacks(x), batchIdx, 'UniformOutput', false);

save([destination filesep 'maskTraceMD.mat'], 'RegStacks', 'batchIdx', 'impath');

% Draw all masks and do it in minibatches
maskTrace = cell(numel(impath), minibatchsize);

% ii is batch number
for ii = 1:numel(impath)
    if  ~all(cellfun(@isempty, maskTrace(ii, :)))
        continue
    end
    
    fprintf('Batch %d/%d', ii, numel(impath));

    imgs = cell(1, numel(impath{ii}));
    
    % Load the img batch
    % jj is image number (resets for each batch)
    for jj = 1:numel(impath{ii})
        if ~isempty(maskTrace{ii, jj})
            continue
        end
        try
            imgs{jj} = squeeze(tiffreadVolume([impath{ii}{jj} filesep dir([impath{ii}{jj} filesep '*norm_reg_Stack.tiff']).name]));
            imgs{jj} = makeViewableImage(imgs{jj}, [8 6 3]); % 755 is ch1-4; this is SHG, FAD, NADH ||  SHG, FAD+HIF, NADH+~HIF
        catch ME
            warning(getReport(ME))
            warning('\n Above error occured while loading from impath{ii}{jj}\n')
        end
    end
    
    % Draw masks
    f1 = figure('Position', [1739 -1455 1280 603]);
    for jj = 1:numel(imgs)
        if ~isempty(maskTrace{ii, jj})
            continue
        end
        fprintf('\n%s', impath{ii}{jj});
        maskTrace{ii, jj} = drawMask(imgs{jj});
        % Save trace
        save([destination filesep 'maskTracesIntermediates.mat'], 'maskTrace', '-mat');
    end
end
close(f1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Conversion of Trace to Mask %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
readyToConvert = ~[cellfun(@isempty, maskTrace)];
for ii = 1:numel(readyToConvert)
    disp(['Batch ' num2str(ii) '/' num2str(sum(readyToConvert(:, 1)))])
    for jj = 1:sum(readyToConvert(ii, :))
        if exist([impath{ii}{jj} filesep 'maskStack.tiff'], 'file')
            continue
        end
        disp(['Image ' num2str(jj) '/' num2str(sum(readyToConvert(ii, :)))])
       try
            iminfo = imfinfo([impath{ii}{jj} filesep 'norm_reg_Stack.tiff']);
            maskSize = [iminfo(1).Height, iminfo(1).Width];
            if isempty(maskTrace{ii, jj})
                continue
            end
            mask = trace2mask(maskTrace{ii, jj}, maskSize);
            TiffSingleWrite(mask, [impath{ii}{jj} filesep 'maskStack.tiff'], 1)
        catch ME
            warning(getReport(ME))
            warning('\n Above error on mask %d from batch %d correspoinding to --> %s \n', jj, ii, impath{ii}{jj})
        end
    end
end
disp('Conversion finished')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Manual selection of mask sign %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maskPath = FolderFinder('maskStack.tiff', 'D:\UMSCC Rad Study');
alreadyDone = FolderFinder('totalMask.tiff', 'D:\UMSCC Rad Study');
maskPath = maskPath(~contains(maskPath, alreadyDone));

% Save one layer masks
for ii = 1:numel(maskPath)
    if numel(imfinfo([maskPath{ii} filesep 'maskStack.tiff'])) == 1
        mask = tiffreadVolume([maskPath{ii} filesep 'maskStack.tiff']);
        imwrite(logical(mask), sprintf('%s%stotalMask.tiff', maskPath{ii}, filesep));
    end
end
alreadyDone = FolderFinder('totalMask.tiff', 'D:\UMSCC Rad Study');
maskPath = maskPath(~contains(maskPath, alreadyDone));

%% Choose mask sign for multilayer masks
% Batch
minibatchsize = 32;
batchIdx = indexBatch(maskPath, minibatchsize, 'shuffle', true);
maskPath = cellfun(@(x) maskPath(x), batchIdx, 'UniformOutput', false);

for ii = 1:numel(maskPath)
    fprintf('Loading batch %d...\n', ii)
    % Load imgs and masks for batch

    imgs = cell(numel(maskPath{ii}), 1);
    masks = cell(numel(maskPath{ii}), 1);
    for jj = 1:numel(maskPath{ii})
        try
            imgs{jj} = makeViewableImage(tiffreadVolume([maskPath{ii}{jj} filesep 'norm_reg_Stack.tiff']), [8 6 3]);
            masks{jj} = tiffreadVolume([maskPath{ii}{jj} filesep 'maskStack.tiff']);
        catch ME
            beep
            imgs{jj} = zeros(1,1);
            masks{jj} = ones(1,1);
            warning(getReport(ME))
        end
    end

    % Choose for batch
    for jj = 1:numel(maskPath{ii})
        fprintf('Mask %d of batch %d: %s\n', jj, ii, maskPath{ii}{jj})
        mask = sum(chooseMask(imgs{jj}, squeeze(masks{jj})), 3);
        imwrite(mask, sprintf('%s%stotalMask.tiff', maskPath{ii}{jj}, filesep));
    end
end

disp('Selection completed')

%%%%%%%%%%%%%%%%%%%
%% AF Processing %%
%%%%%%%%%%%%%%%%%%%
RegStacks = FolderFinder('norm_reg_Stack.tiff', [destination filesep '..' filesep 'Stitched']);
AF = RegStacks([contains(RegStacks, 'AF')]);
np = pathpieces(AF);
orrTab = cell(1+numel(AF), 12);
orrTab(1, :) = {'Tx', 't', 'Animal', 'Sample', 'average pixel-wise ORR', 'average FAD Int', 'average NAD(P)H int', 'bulk ORR', 'var pixel-wise ORR ', 'var average FAD Int', 'var average NAD(P)H int', 'n'};
orrTab(2:end, 1:4) = np(:, 5:8);
% errors = cell(0, 2);

for ii = 1:numel(AF)
    sd = strrep(AF{ii}, 'Stitched','Processed');
%     if exist(sprintf('%s\\orrMap.tiff', sd), 'file') && exist(sprintf('%s\\orrMapMasked.tiff', sd), 'file') && exist(sprintf('%s\\avgIntensity.tiff', sd), 'file') && exist(sprintf('%s\\orrFig.png', sd), 'file') && exist(sprintf('%s\\orrFigMasked.png', sd), 'file')
%         fprintf('Skipping %s', AF{ii})
%         continue
%     end
    try
        fprintf('Pair %d/%d >>> %s\n\n', ii, numel(AF), AF{ii});
                
        % Load Image
        img = squeeze(tiffreadVolume([AF{ii} filesep 'norm_reg_Stack.tiff']));
        mask = logical(squeeze(tiffreadVolume([AF{ii} filesep 'totalMask.tiff'])));

        % Check if the mask is just a placeholder for a bad mask
        if all(mask, 'all') || all(~mask, 'all')
            warning('Suspected bad img based on empty mask for %s. Manually verify and rerun if necessary.\n', AF{ii})
            continue
        end

        % Calculate Redox
        nadh = medfilt2(img(:,:,3), [5 5]);  % 755 Blue
        fad = medfilt2(img(:,:,6), [5 5]);   % 855 Green
        % Add thresholds to the mask, get rid of background
        mask = mask & fad > 3e-4 & nadh > 3e-4;
        int = mask.*(fad + nadh)/2;
        orr = mask.*(fad./(fad+nadh));
        cORR = createColorORRMap(orr, int, 'redoxMax', 0.8);

        % Rewrite Updated mask
        imwrite(mask, [AF{ii} filesep 'totalMask.tiff']);

        % Save outputs to table (these are all masked)
        r = ii+1;
        orrTab(r, 5:7) = cellfun(@(x) mean(x, 'all', 'omitnan'), {orr(mask), fad(mask), nadh(mask)}, 'UniformOutput', false);
        orrTab{r, 8} = orrTab{r,6}/(orrTab{r, 6} + orrTab{r,7});
        orrTab(r, 9:11) = cellfun(@(x) var(x, 0, 'all', 'omitnan'), {orr(mask), fad(mask), nadh(mask)}, 'UniformOutput', false);
        orrTab{r, 12} = sum(mask, 'all');

        
        % Save images/figures
        if ~exist(sd, 'dir')
            mkdir(sd);
        end
        TiffSingleWrite(orr, sprintf('%s\\orrMap.tiff', sd), 1)
        TiffSingleWrite(orr.*mask, sprintf('%s\\orrMapMasked.tiff', sd), 1);
        TiffSingleWrite(int, sprintf('%s\\avgIntensity.tiff', sd), 1);
        imwrite(cORR, sprintf('%s\\orrFig.png', sd), 'png');
        imwrite(cORR.*uint8(mask), sprintf('%s\\orrFigMasked.png', sd), 'png');

    catch ME
        beep
        errors{end+1, 1} = AF{ii};
        errors{end, 2} = getReport(ME);
        warning('\n %s \n', getReport(ME))
        fprintf('\nAbove error on %s\n\n', AF{ii})
    end
end

writecell(orrTab, [destination filesep 'orrTable.xlsx'])
beep
disp('ORR Calcs complete')

%%%%%%%%%%%%%%%%%%%%
%% IHC Processing %%
%%%%%%%%%%%%%%%%%%%%
%%% notes on thresholds %%%
% LP1=2.2859 is the graythresh value for all the ihc pixels in the set of
% IHC images with AF pairing.
% The quality of these thressholds was checked on multiple random samples.
% LH=10 was experimentally determined using graythresh values and histogram
% analysis as a guide

IHC  = RegStacks([contains(RegStacks, 'IHC')]);
np = pathpieces(IHC);
ihcTab = cell(1+numel(IHC), 11);
ihcTab(1, :) = {'Tx', 't', 'Animal', 'Sample', 'hif regional average', 'pimo regional average', 'total hif fraction', 'total pimo fraction', 'hif regional var', 'pimo regional var', 'n'};
ihcTab(2:end, 1:4) = np(:, 5:8);

% Set threshholds
H = fspecial('average', 512);
LP = 2.2859;
LH = 10;

% Filter prep
fR = 69 + (1/3);   % This will give us a 100 um disk (the general O2 limit)
fA = pi*fR^2;
f = fA*fspecial('disk', fR);

for ii = 1:numel(IHC)
    sd = strrep(IHC{ii}, 'Stitched','Processed');
%     if exist(sprintf('%s\\pim_hif_shg_Map.tiff', sd), 'file') && exist(sprintf('%s\\PimovHifMaskedFilteredMasked.tiff',sd), 'file') && exist(sprintf('%s\\PimovHifImage.tiff', sd), 'file') && exist(sprintf('%s\\PimovHifFig.tiff', sd), 'file')
%         fprintf('Skipping %s', IHC{ii})
%         continue
%     end
    try
        fprintf('\n Pair %d/%d >>> %s\n\n', ii, numel(IHC), IHC{ii});
    
        % Load Image
        img = squeeze(tiffreadVolume([IHC{ii} filesep 'norm_reg_stack.tiff']));
        mask = logical(tiffreadVolume([IHC{ii} filesep 'totalMask.tiff']));

        % Check if the mask is just a placeholder for a bad mask
        if all(mask, 'all') || all(~mask, 'all')
            warning('Suspected bad img based on empty mask for %s. Manually verify and rerun if necessary.\n', IHC{ii})
            continue
        end

        % Extract signals
        pim = img(:,:,1);   % 755 Red
        pimPos = pim>LP;
        hif = img(:,:,5).*img(:,:,6); % 950 Red * 950 Green
        hifPos = hif>LH;
        shg = img(:,:,7); % 950 Blue

        % Regional sampling
        PvH = mask.*imfilter(double(mask.*cat(3, pimPos, hifPos)), f);
        reg = imfilter(double(mask), f);
        PvH = PvH./reg; % Convert pixel counts to fraction of pixels in the region
        pimReg = PvH(:,:,1);
        hifReg = PvH(:,:,2);
        PvHIm = cat(3, PvH, zeros(size(PvH, [1 2])));
        
        % Save Outputs to table
        r = ii+1;
        ihcTab(r, 5:6) = cellfun(@(x) mean(x, 'all', 'omitnan'), {hifReg(mask), pimReg(mask)}, 'UniformOutput', false);
        ihcTab(r, 7:8) = cellfun(@(x) sum(x, 'all', 'omitnan')/numel(x), {hifPos(mask), pimPos(mask)}, 'UniformOutput', false);
        ihcTab(r, 9:10) = cellfun(@(x) var(x, 0, 'all', 'omitnan'), {hifReg(mask), pimReg(mask)}, 'UniformOutput', false);
        ihcTab{r, 11} = sum(mask, 'all');

         % Save images/figures

        if ~exist(sd, 'dir')
            mkdir(sd);
        end
        TiffSingleWrite(cat(3, pim, hif, shg), sprintf('%s\\pim_hif_shg_Map.tiff', sd), 1);
        TiffSingleWrite(PvH, sprintf('%s\\PimovHifMaskedFilteredMasked.tiff',sd), 1);
        TiffSingleWrite(PvHIm, sprintf('%s\\PimovHifImage.tiff', sd), 1);

    catch ME
        beep
        errors{end+1, 1} = IHC{ii};
        errors{end, 2} = getReport(ME);
        warning(getReport(ME))
        fprintf('\nAbove error on %s\n\n', IHC{ii})
    end
end
writecell(ihcTab, [destination filesep 'ihcTable.xlsx'])
beep
disp('IHC Finished')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Intermodal Manual Registration %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check which samples have both modes (for multimodal analysis);
% If bad image(s) appear in pair: exit out of image figure. This will throw
% and error that we catch and will display a warning confirming it. It will
% also add the images to a list for "manualRework".

pp = pathpieces(RegStacks);
bothProc = {};
for ii= 1:size(pp, 1)
    sam = namemaker(pp(ii,1:end-2), filesep);
    if mod(numel(RegStacks(contains(RegStacks, sam))), 2) == 0 ...
        && sum(contains(RegStacks(contains(RegStacks, sam)), 'AF')) == 0.5*numel(RegStacks(contains(RegStacks, sam))) ...
        && sum(contains(RegStacks(contains(RegStacks, sam)), 'IHC')) == 0.5*numel(RegStacks(contains(RegStacks, sam)))
        bothProc{end+1} = RegStacks(contains(RegStacks, sam));
    end
end
staq = vertcat(bothProc{:});
bothProc = unique(staq(:,1));
bothProc(:,2) = unique(staq(:,2));

% Batch
minibatchsize = 16;
batchIdx = indexBatch(bothProc(:,1), minibatchsize, 'shuffle', true);
bothProc = cellfun(@(x) bothProc(x, :), batchIdx, 'UniformOutput', false);

transforms = cell(numel(bothProc), 1);

% Find saved normalized image pairs
for ii = 1:numel(bothProc)
        % Load the batch
        orr = cell(size(bothProc{ii}, 1), 1);
        phs = cell(size(bothProc{ii}, 1), 1);
        transforms{ii} = cell(size(bothProc{ii}, 1), 1);

        % Load the images (if they haven't already been registered)
        for jj = 1:size(bothProc{ii}, 1)
            try
                if isempty(transforms{ii}{jj}) && ~exist(sprintf('%s\\..\\..\\reg_trans.mat', bothProc{ii}{jj, 1}), 'file')
                    % Load image pair
                    orr{jj} = squeeze(tiffreadVolume([bothProc{ii}{jj, [contains(bothProc{ii}(jj,:), 'AF')]} filesep 'norm_reg_Stack.tiff']));
                    phs{jj} = squeeze(tiffreadVolume([bothProc{ii}{jj, [contains(bothProc{ii}(jj,:), 'IHC')]} filesep 'norm_reg_Stack.tiff']));
                else
                    fprintf('reg_trans.mat previously created for %s...skipped\n', bothProc{ii}{jj, 1})
                end
            catch ME
                beep
                orr{jj} = zeros(1);
                phs{jj} = zeros(1);
                errors{end+1, 1} = bothProc{ii}{jj, 1};
                errors{end, 2} = getReport(ME);
                warning('\n %s', getReport(ME))
                fprintf('\nAbove error on %s\n\n', bothProc{ii}{jj, 1})
            end
        end
        
        % Register the images that were loaded
        for jj = 1:numel(orr)
            fprintf('Currently registering %s\n', bothProc{ii}{jj});
            try
                % Get registration transform (ihc is moving: ch5-8)
                [~, transforms{ii}{jj}] = alignByLineROI(makeViewableImage(orr{jj}, [4 2 3]), makeViewableImage(phs{jj}, [7 5 1]));
                t = transforms{ii}{jj};
                save(sprintf('%s\\..\\..\\reg_trans.mat', bothProc{ii}{jj, 1}), 't', '-mat')
            catch ME
                switch ME.identifier
                    case 'MATLAB:hg;DeletedObject'
                        manualRework{end+1} = bothProc{ii}{jj};
                    case 'MATLAB:badsubscript'
                        warning('No image loaded for %s', bothProc{ii}{jj});
                    otherwise
                        errors{end+1, 1} = bothProc{ii}{jj, 1};
                        errors{end, 2} = getReport(ME);
                        warning('\n %s', getReport(ME))
                        fprintf('\nAbove error on %s\n\n', bothProc{ii}{jj, 1})
                end
            end
        end
end
bothProc = vertcat(bothProc{:}); % Unbatch
disp('Complete.')
rt = FolderFinder('reg_trans.mat');
%% Check for transforms that were not saved properly
for ii = 1:numel(rt)
    mf = matfile(sprintf('%s\\reg_trans.mat', rt{ii}));
    if isempty(mf.t)
        try
            af = FolderFinder('norm_reg_Stack.tiff', [rt{ii} filesep 'AF']);
            ihc = FolderFinder('norm_reg_Stack.tiff', [rt{ii} filesep 'IHC']);
            orr = squeeze(tiffreadVolume([af{:} filesep 'norm_reg_Stack.tiff']));
            phs = squeeze(tiffreadVolume([ihc{:} filesep 'norm_reg_Stack.tiff']));
        catch ME
            warning('%s \n \n %s failed to load\n', getReport(ME), rt{ii})
        end
         % Get registration transform (ihc is moving: ch5-8)
         try
            [~, t] = alignByLineROI(makeViewableImage(orr, [4 2 3]), makeViewableImage(phs, [7 5 1]));
            save(sprintf('%s\\reg_trans.mat', rt{ii}), 't', '-mat')
         catch ME
         end
    else
        disp([rt{ii} 'not empty'])
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Intermodal Processing %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
np = pathpieces(bothProc(:,1));
imTab = cell(1+size(bothProc, 1), 11);
imTab(1, :) = {'Tx', 't', 'Animal', 'Sample',  'avg pixel-wise ORR', 'avg pimo regional', 'avg hif regional', 'var pixel-wise ORR', 'var pimo regional', 'var hif regional', 'total region pixel count'};
imTab(2:end, 1:4) = np(:, 5:8);

% Filter prep
fR = 69 + (1/3);   % This will give us a 100 um disk (the O2 limit)
fA = pi*fR^2;
f = fA*fspecial('disk', fR); % Each pixel will count as one (except edges which are weighted less than 1), then divideing by the region size should give the percentage at each region
for ii = 1:size(bothProc, 1)
    try
        load(sprintf('%s\\..\\..\\reg_trans.mat', bothProc{ii, 1}), 't');
        orr = squeeze(tiffreadVolume(strrep([bothProc{ii, contains(bothProc(ii, :), 'AF')} filesep 'orrMap.tiff'], 'Stitched', 'Processed')));
        mask = logical(squeeze(tiffreadVolume(([bothProc{ii, contains(bothProc(ii, :), 'AF')} filesep 'totalMask.tiff']))));
        phs = squeeze(tiffreadVolume(strrep([bothProc{ii, contains(bothProc(ii, :), 'IHC')} filesep 'pim_hif_shg_Map.tiff'], 'Stitched', 'Processed')));
        ihcMask = logical(squeeze(tiffreadVolume(([bothProc{ii, contains(bothProc(ii, :), 'IHC')} filesep 'totalMask.tiff']))));

        % Transform outputs by registration guidelines
        phs = imageTransformer(phs, t);
        ihcMask = imageTransformer(ihcMask, t);

        % Pad to match size and properly overlap
        % Images
        size2pad = max(cell2mat(cellfun(@(x) size(x, [1 2 3]), {orr; phs}, 'UniformOutput', false))); % Find largest dims overall
        orr(end+1:size2pad(1), end+1:size2pad(2)) = 0;
        phs(end+1:size2pad(1), end+1:size2pad(2), :) = 0;
        pairOut = cat(3, orr, phs);

        % Masks
        size2pad = max(cell2mat(cellfun(@size, {mask; ihcMask}, 'UniformOutput', false))); % Find largest dims overall
        mask(end+1:size2pad(1), end+1:size2pad(2)) = 0;
        ihcMask(end+1:size2pad(1), end+1:size2pad(2)) = 0;
        pairMask = mask & ihcMask;

        % Signal extraction
        pim = phs(:,:,1);
        pimPos = pim>LP;
        hif = phs(:,:,2);
        hifPos = hif>LH;
        sigs = cat(3, orr, pimPos, hifPos);
        
        % Intermodal regional analysis
        ophOut = pairMask.*imfilter(double(pairMask.*sigs), f);
        regionSize = pairMask.*imfilter(double(pairMask), f);
        regAvg = ophOut./regionSize;

        % Save outputs to table
        r = ii+1;
        imTab(r, 5:7) = num2cell(mean(regAvg, [1 2], 'omitnan'));
        imTab(r, 8:10) = num2cell(var(regAvg, 0, [1 2], 'omitnan'));
        imTab{r, 11} = sum(regionSize, 'all')/fA;

        % Save images
        TiffSingleWrite(pairOut, [bothProc{ii, 1} filesep '..' filesep '..' filesep 'intermodalRegisteredStack.tiff'], 1)
        TiffSingleWrite(pairMask, [bothProc{ii, 1} filesep '..' filesep '..' filesep 'intermodalRegisteredMask.tiff'], 1)
        TiffSingleWrite(regAvg, [bothProc{ii, 1} filesep '..' filesep '..' filesep 'orr_pimo_hif_RegionalAnalysis.tiff'], 1)

        fprintf('>>> %d Complete\n', ii)
    catch ME
        beep
        errors{end+1, 1} = bothProc{ii, 1};
        errors{end, 2} = getReport(ME);
        warning('\n %s', getReport(ME))
        fprintf('\nAbove error on %s\n\n', bothProc{ii, 1})
    end
end

writecell(imTab, [destination filesep 'intermodalTable.xlsx'])
disp('Intermodal processing complete');

warning('on', 'MATLAB:imagesci:Tiff:missingExtraSamples')

%% Make intermodal figures
IMProc = FolderFinder('orr_pimo_hif_RegionalAnalysis.tiff', [destination filesep '..' filesep 'Stitched']);
pp = pathpieces(IMProc);
cat1 = unique(pp(:, 4));
cat2 = unique(pp(:, 5));
cat3 = unique(pp(:, 6));

%%
% Create and save cummulative phasors
fig = figure;
fig.WindowState = 'fullscreen';
ax = axes(fig);
ax.FontName = 'arial';
ax.FontSize = 18;
cmap = colormap('jet').*colormap('gray');

for ii = 1:numel(cat1)
    sub1 = [destination filesep cat1{ii}];
    
    % Pimo 1d Hist Fig & ax
    phis = figure('WindowState', 'fullscreen');
    clf(phis)
    phax = axes(phis);
    hold(phax, 'on');

    % HIF 1d Hist Fig & ax
    hhis = figure('WindowState', 'fullscreen');
    clf(hhis)
    hhax = axes(hhis);
    hold(hhax, 'on');

    % ORR 1d Hist Fig & ax
    ohis = figure('WindowState', 'fullscreen');
    clf(ohis)
    ohax = axes(ohis);
    hold(ohax, 'on');

    % Linestyles
    c = cmap([0.25 0.5 0.75]*256, :);
    l = {'-', ':', '-.'};
    m = {'square', 'o', 'v'};
    ls = 0;
    

    for jj = 1:numel(cat2)
        sub2 = [sub1 filesep cat2{jj}];
        for kk = 1:numel(cat3)
            if ~(strcmp(cat2{jj}, 'XT') && (strcmp(cat3{kk}, '24hpt') || strcmp(cat3{kk}, '48hpt'))) && ~strcmp(cat3{kk}, '0hpt-Baseline')
                continue
            end
            sub3 = [sub2 filesep cat3{kk}];
            ims = cellfun(@(x) squeeze(tiffreadVolume([x filesep 'orr_pimo_hif_RegionalAnalysis.tiff'])), IMProc(strcmp(pp(:,4), cat1{ii}) & strcmp(pp(:,5), cat2{jj}) & strcmp(pp(:,6), cat3{kk})), 'UniformOutput', false);
            if isempty(ims)
                continue
            end
            fprintf('%s %s %s n=%d \n', cat1{ii}, cat2{jj}, cat3{kk}, numel(ims))
            ims = cellfun(@(x) reshape(x, prod(size(x, [1 2])), []), ims, 'UniformOutput', false);
            ims = vertcat(ims{:});
            bad = any(isnan(ims), 2) | all(ims==0, 2);
            orr = ims(~bad,1);
            pim = 100*ims(~bad,2); % Convert to %
            hif = 100*ims(~bad,3); % Convert to %

            % Prep bins for each
            nbins = 25;
            pimLim = [0 100];
            hifLim = [0 100];
            orrLim = [0 1];

            P = linspace(pimLim(1), pimLim(2), nbins+1);
            H = linspace(hifLim(1), hifLim(2), nbins+1);
            O = linspace(orrLim(1), orrLim(2), nbins+1);

            if (strcmp(cat2{jj}, 'XT') && (strcmp(cat3{kk}, '24hpt') || strcmp(cat3{kk}, '48hpt'))) || strcmp(cat3{kk}, '0hpt-Baseline')
                fprintf('Histing %s-%s-%s', cat1{ii}, cat2{jj}, cat3{kk})
                ls = ls+1;
                spec = {'LineWidth', 5, 'Color', c(ls, :), 'Marker', m{ls}, 'Linestyle', l{ls}};

                % 1D Hist of ORR
                N = histcounts(orr(orr>=orrLim(1) & orr<=orrLim(2)), O);
                plot(ohax, O(1:end-1)+(diff(orrLim)/nbins/2), 100*(N/numel(orr)), spec{:})
    
                ohax.FontName = 'arial';
                ohax.FontSize = 30;
                ohax.XLabel.String = 'Mean ORR';
                ohax.YLabel.String = '% of Pixels';
                ohax.YScale = 'linear';
                ohax.YLim = [0 35];
                ohax.YTick = 0:10:30;

                % 1D Hist of Pimo
                N = histcounts(pim(pim>=pimLim(1) & pim<=pimLim(2)), P);
                plot(phax, P(1:end-1)+(diff(pimLim)/nbins/2), 100*(N/numel(pim)), spec{:})
    
                phax.FontName = 'arial';
                phax.FontSize = 30;
                phax.XLabel.String = 'Pimo+%';
                phax.YLabel.String = '% of Pixels';
                phax.YScale = 'log';
                phax.YLim = [10^-3 10^2];
    
                % 1D Hist of HIF
                N = histcounts(hif(hif>=hifLim(1) & hif<=hifLim(2)), H);
                plot(hhax, H(1:end-1)+(diff(hifLim)/nbins/2), 100*(N/numel(hif)), spec{:})
                hhax.FontName = 'arial';
                hhax.FontSize = 30;
                hhax.XLabel.String = 'HIF1\alpha+%';
                hhax.YLabel.String = '% of Pixels';
                hhax.YScale = 'log';
                hhax.YLim = [10^-3 10^2];
            end

            P = linspace(pimLim(1)+0.5*diff(pimLim)/nbins, pimLim(2)-0.5*diff(pimLim)/nbins, nbins);
            H = linspace(hifLim(1)+0.5*diff(hifLim)/nbins, hifLim(2)-0.5*diff(hifLim)/nbins, nbins);
            O = linspace(orrLim(1)+0.5*diff(orrLim)/nbins, orrLim(2)-0.5*diff(orrLim)/nbins, nbins);

            % Pimo v HIF
            ctrs = {P, H};
            data = [pim(pim>pimLim(1) & pim<=pimLim(2) & hif>hifLim(1) & hif<=hifLim(2)), hif(pim>pimLim(1) & pim<=pimLim(2) & hif>hifLim(1) & hif<=hifLim(2))];
            
            % Plot
            [N, ctrs] = hist3(data, 'Ctrs', ctrs);
            N = log10(100*N'/size(data, 1));
            surf(ax, ctrs{:}, N);
            
            shading(ax, 'interp')
            view(ax, 2)
            ax.Title.String = sprintf('%s %s %s', cat1{ii}, cat2{jj}, cat3{kk});
            ax.Title.FontSize = 28;
            ax.FontName = 'arial';
            ax.FontSize = 18;
            ax.XLabel.String = 'Pimo+%';
            ax.YLabel.String = 'HIF1\alpha+%';
            ax.XLabel.FontSize = 24;
            ax.YLabel.FontSize = 24;
            ax.XLim = [min(ctrs{1}), max(ctrs{1})];
            ax.YLim = [min(ctrs{2}), max(ctrs{2})];
            colormap(ax, cmap)
            clim(ax, [-5 1])
            cb = colorbar(ax);
            cb.Label.String = 'Log_{10} % of Pixels';
            cb.Label.FontSize = 24;
            savefig(fig, sprintf('%s%spimo-%.2f-%.2f_hif-%.2f-%.2f_2dHist.fig', sub3, filesep, pimLim, hifLim));
            exportgraphics(fig, sprintf('%s%spimo-%.2f-%.2f_hif-%.2f-%.2f_2dHist.png', sub3, filesep, pimLim, hifLim), 'Resolution', 300);
            
            % Pimo v ORR
            ctrs = {P, O};
            data = [pim(pim>pimLim(1) & pim<=pimLim(2) & orr>orrLim(1) & orr<1), orr(pim>pimLim(1) & pim<=pimLim(2) & orr>orrLim(1)& orr<1)];
            
            % Plot
            [N, ctrs] = hist3(data, 'Ctrs', ctrs);
            surf(ax, ctrs{:}, log10(100*N'/size(data, 1)));
            
            shading(ax, 'interp')
            view(ax, 2)
            ax.Title.String = sprintf('%s %s %s', cat1{ii}, cat2{jj}, cat3{kk});
            ax.Title.FontSize = 28;
            ax.FontName = 'arial';
            ax.FontSize = 18;
            ax.XLabel.String = 'Pimo+%';
            ax.YLabel.String = 'Mean ORR';
            ax.XLabel.FontSize = 24;
            ax.YLabel.FontSize = 24;
            ax.XLim = [min(ctrs{1}), max(ctrs{1})];
            ax.YLim = [min(ctrs{2}), max(ctrs{2})];
            colormap(ax, cmap)
            clim(ax, [-5 1])
            cb = colorbar(ax);
            cb.Label.String = 'Log_{10} % of Pixels';
            cb.Label.FontSize = 24;
            savefig(fig, sprintf('%s%spimo-%.2f-%.2f_orr-%.2f-%.2f_2dHist.fig', sub3, filesep, pimLim, orrLim));
            exportgraphics(fig, sprintf('%s%spimo-%.2f-%.2f_orr-%.2f-%.2f_2dHist.png', sub3, filesep, pimLim, orrLim), 'Resolution', 300);

            % HIF v ORR
            ctrs = {H, O};
            data = [hif(hif>hifLim(1) & hif<=hifLim(2) & orr>orrLim(1) & orr<=orrLim(2)), orr(hif>hifLim(1) & hif<=hifLim(2) & orr>orrLim(1) & orr<=orrLim(2))];
            
            % Plot
            [N, ctrs] = hist3(data, 'Ctrs', ctrs);
            surf(ax, ctrs{:}, log10(100*N'/size(data, 1)));
            
            shading(ax, 'interp')
            view(ax, 2)
            ax.Title.String = sprintf('%s %s %s', cat1{ii}, cat2{jj}, cat3{kk});
            ax.Title.FontSize = 28;
            ax.FontName = 'arial';
            ax.FontSize = 18;
            ax.XLabel.String = 'HIF1\alpha+%';
            ax.YLabel.String = 'Mean ORR';
            ax.XLabel.FontSize = 24;
            ax.YLabel.FontSize = 24;
            ax.XLim = [min(ctrs{1}), max(ctrs{1})];
            ax.YLim = [min(ctrs{2}), max(ctrs{2})];
            colormap(ax, cmap)
            clim(ax, [-5 1])
            cb = colorbar(ax);
            cb.Label.String = 'Log_{10} % of Pixels';
            cb.Label.FontSize = 24;
            savefig(fig, sprintf('%s%shif-%.2f-%.2f_orr-%.2f-%.2f_2dHist.fig', sub3, filesep, hifLim, orrLim));
            exportgraphics(fig, sprintf('%s%shif-%.2f-%.2f_orr-%.2f-%.2f_2dHist.png', sub3, filesep, hifLim, orrLim), 'Resolution', 300);

            % Trivariate Histogram
            ctrs = {P, H};
            data = [pim(pim>pimLim(1) & pim<=pimLim(2) & hif>hifLim(1) & hif<=hifLim(2) & orr>orrLim(1) & orr<=orrLim(2)), ...
                    hif(pim>pimLim(1) & pim<=pimLim(2) & hif>hifLim(1) & hif<=hifLim(2) & orr>orrLim(1) & orr<=orrLim(2)), ...
                    orr(pim>pimLim(1) & pim<=pimLim(2) & hif>hifLim(1) & hif<=hifLim(2) & orr>orrLim(1) & orr<=orrLim(2))];
            N = hist3(data(:,1:2), 'Ctrs', ctrs);
            C = HistogramColorCalculator(data(:,1), data(:,2), data(:,3),  ctrs);
            surf(ax, ctrs{:}, log10(100*N'/size(data, 1)), C);
            shading(ax, 'interp')
            view(ax, 45, 60)
            ax.XLabel.String = 'Pimo+%';
            ax.YLabel.String = 'HIF1\alpha+%';
            ax.ZLabel.String = 'Log_{10} % of Pixels';
            ax.XLabel.FontSize = 24;
            ax.YLabel.FontSize = 24;
            ax.ZLabel.FontSize = 24;
            colormap(ax, cmap)
            clim(ax, [0 0.8])
            cb = colorbar(ax);
            cb.Label.String = 'Mean ORR';
            cb.Label.FontSize = 24;
            savefig(fig, sprintf('%s%spimo-%.2f-%.2f_hif-%.2f-%.2f_3dHistColored.fig', sub3, filesep, hifLim, pimLim));
            exportgraphics(fig, sprintf('%s%spimo-%.2f-%.2f_hif-%.2f-%.2f_3dHistColored.png', sub3, filesep, hifLim, pimLim), 'Resolution', 300);
%             
        end
    end
    legend(ohax, 'Baseline', '24hpt', '48hpt', 'Location', 'northeastoutside')
    legend(phax, 'Baseline', '24hpt', '48hpt', 'Location', 'northeastoutside');
    legend(hhax, 'Baseline', '24hpt', '48hpt', 'Location', 'northeastoutside');
    savefig(ohis, sprintf('%s\\%s\\ORR_1DHistogram-LinearScale-%.2f-%.2f.fig', destination, cat1{ii}, pimLim));
    exportgraphics(ohax, sprintf('%s\\%s\\ORR_1DHistogram-LinearScale-%.2f-%.2f.png', destination, cat1{ii}, pimLim), 'Resolution', 300);
    savefig(phis, sprintf('%s\\%s\\PimoPos_1DHistogram-LinearScale-%.2f-%.2f.fig', destination, cat1{ii}, pimLim))
    exportgraphics(phax, sprintf('%s\\%s\\PimoPos_1DHistogram-LinearScale-%.2f-%.2f.png', destination, cat1{ii}, pimLim), 'Resolution', 300);
    savefig(hhis, sprintf('%s\\%s\\HIFPos_1DHistogram-LinearScale-%.2f-%.2f.fig', destination, cat1{ii}, hifLim))
    exportgraphics(hhax, sprintf('%s\\%s\\HIFPos_1DHistogram-LinearScale-%.2f-%.2f.png', destination, cat1{ii}, hifLim), 'Resolution', 300);
end

%% Group and individual stats
repPaths = {};
for ii = 1:numel(cat1)
    sub1 = [destination filesep cat1{ii}];
    for jj = 1:numel(cat2)
        sub2 = [sub1 filesep cat2{jj}];
        for kk = 1:numel(cat3)
            sub3 = [sub2 filesep cat3{kk}];
            groupIms = IMProc(strcmp(pp(:,4), cat1{ii}) & strcmp(pp(:,5), cat2{jj}) & strcmp(pp(:,6), cat3{kk}));
            if isempty(groupIms)
                continue
            end
            ims = cellfun(@(x) squeeze(tiffreadVolume([x filesep 'orr_pimo_hif_RegionalAnalysis.tiff'])), groupIms, 'UniformOutput', false);
            individualStats = cellfun(@(x) [squeeze(mean(x, [1 2], 'omitnan')); squeeze(var(x, 1, [1 2], 'omitnan'))], ims, 'UniformOutput', false);

            ims = cellfun(@(x) reshape(x, [prod(size(x, [1 2])), 3]), ims, 'UniformOutput', false);
            ims = vertcat(ims{:});
            overallStat = [squeeze(mean(ims, 1, 'omitnan')), squeeze(var(ims, 1, 1, 'omitnan'))];

            rep = cellfun(@(x) abs(x'-overallStat), individualStats, 'UniformOutput', false);
            rep = mean(vertcat(rep{:}), 2);
            repPaths(end+1) = groupIms(rep==min(rep));
            fprintf('%s\n', repPaths{end});
        end
    end
end