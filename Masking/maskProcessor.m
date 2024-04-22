%% -- Mask Processor -- %%
% % Flow of Program % %
% This program is divided into parts to optimize the efficiency of manual
% mask creation. First, the images will be displayed and you will be
% allowed to trace masks on each image. This is done in batches, so loading
% time for a single batch is long, but the time between images within the
% batch is minimal. This allows for mask-tracing in "waves." 
% 
% Next, the masks traces will be converted into masks (logical polygons).
% This process is computationally intensive, ut entirely input free,
% meaning you, the user, can be absent during this process. All ouptuts
% from this stage will be written to the HD, so the next step can be
% returned to even after MATLAB has closed. You will recieve a message in
% the command prompt when this conversion is finished.
%
% Finally, masks are selected. This is when you will shoose whether or not
% each trace is a positive or negative mask. This processi s batched
% similarly to part 1, allowing selection in waves with long gaps for
% laoding. Each batch is written to the HD, and then the next batch is
% loaded.


% % Instructions for use % %
% In the tracing stage, you will click to begin a trace, then duble click
% to complete it. Click again to prepared for another trace. Click to start
% it and double click to complete it. This can be repeated as many times as
% necessary. When you ahve tracd all that is needed, press enter to move to
% the next image.
%
% No input from the user is required for the conversion stage.
%
% In the selection stage, the image with the first part (from the first
% trace) of the mask will be presented. Press space to flip the mask or
% enter to accept it and move to the next section of the mask. In many
% cases, when moving to the next seciton, nothing about the displaywill
% change UNTIL YOU HIT SPACE, so make sure once you have hit enter, you hit
% space to check which part of the mask is under consdieration. You can hit
% space multiple times to alternate the mask and compare choices. 

% % Special Marker Keys During Selection Stage % %
% In the selection phase, you have the option to mark masks/images for
% manual evaluation using the following keys:
% F -- Image is noisy (but otherwise acceptable). Pressing F will not
% impact your ability to choose mask sign. You will still choose as normal,
% but a note will be saved with the mask.
%
% X/Esc -- Image is unacceptable (bad alignment of channels, bad stitching,
% etc.) Pressing X/Esc will mark the image and advance you to the next
% image.

%% Mask tracing
% Where are we looking?
destination = 'D:\UMSCC Rad Study\Processed';

% Get registered stack locations
RegStacks = FolderFinder('intermodalRegisteredStack', destination);

% minibatching to improve time between sets of images. Shuffling to limit
% any systematic error in tracing from working at different times
minibatchsize = 44;
batchIdx = indexBatch(RegStacks, minibatchsize, 'shuffle', true);
impath = cellfun(@(x) RegStacks(x), batchIdx, 'UniformOutput', false);
%%
save([destination filesep 'maskTraceMultiModeMD.mat'], 'RegStacks', 'batchIdx', 'impath');

% Draw all masks and do it in minibatches
maskTrace = cell(numel(impath), minibatchsize);
for ii = 1:numel(impath)
    imgs = cell(1, numel(impath{ii}));
    
    % Load the img batch
    for jj = 1:numel(impath{ii})
        imgs{jj} = squeeze(tiffreadVolume([impath{ii}{jj} filesep dir([impath{ii}{jj} filesep '*registeredStack.tiff']).name]));
        if size(imgs{jj}, 3) == 8
            imgs{jj} = makeViewableImage(imgs{jj}, [8 6 3]);
        else
            imgs{jj} = makeViewableImage(imgs{jj}, [4 2 6]);
        end
    end
    
    % Draw masks
    maskTrace{ii} = drawMask(imgs);
end

% Save intermediate
save(sprintf('%s%smaskTraceMultiMode.mat', destination, filesep), 'maskTrace');

%% Convert all masks
maskTrace = maskTrace(:,1);
noTrace = cellfun(@(x) cellfun(@isempty, x), maskTrace, 'UniformOutput', false);
emptMask = {};
errors = cell(0, 2);
for ii = 1:numel(impath)
    try
    disp(['Batch ' num2str(ii) '/' num2str(numel(impath))])
    for jj = 1:sum(cellfun(@(x) ~isempty(x), impath{ii}))
        disp(['Image ' num2str(jj) '/' num2str(sum(cellfun(@(x) ~isempty(x), impath{ii})))])
        if exist([impath{ii}{jj} filesep 'maskStack.tiff'], 'file')
            continue
        end
        iminfo = imfinfo([impath{ii}{jj} filesep dir([impath{ii}{jj} filesep '*egisteredStack.tiff']).name]);
        maskSize = [iminfo(1).Height, iminfo(1).Width];
        trace = maskTrace{ii}(jj, ~noTrace{ii}(jj,:));
        if isempty(trace)
            emptMask{end+1} = impath{ii}{jj};
            continue
        end
        mask = trace2mask(trace, maskSize);
        TiffSingleWrite(mask, [impath{ii}{jj} filesep 'maskStack.tiff'], 1)
    end
    catch ME
        warning(getReport(ME))
        errors{end+1, 1} = impath{ii}{jj};
        errors{end, 2} = getReport(ME);
    end
end

disp('Conversion finished')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you want to stop session after drawing a mask, use "Run Section"     %
% instead of "Run". The intermediate masks (before selection of sign) are %
% saved and will be reaccessed in the following section                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: Add something to deal with images masks are marked                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Choose mask sign
% Find and batch mask stacks
maskPath = FolderFinder('maskStack.tiff', destination);
maskSeld = FolderFinder('totalMask.tiff', destination);
maskPath = maskPath(~cellfun(@(x) any(strcmp(maskSeld, x)), maskPath));
minibatchsize = 44;
batchIdx = indexBatch(maskPath, minibatchsize, 'shuffle', true);
impath = cellfun(@(x) maskPath(x), batchIdx, 'UniformOutput', false);
notes = cell(0, 3); % Cell array for noting problematic and/or noisy images for manual evaluation

%% Batch through and choose
for ii = 1:numel(impath)
    fprintf('\n Batch %d/%d', ii, numel(impath));
    imgs = cell(numel(impath{ii}), 1);
    masks = cell(numel(impath{ii}), 1);
    mask = cell(numel(masks), 1);

    try
        % Load batch
        for jj = 1:numel(impath{ii})
            imgs{jj} = makeViewableImage(tiffreadVolume([impath{ii}{jj} filesep dir([impath{ii}{jj} filesep '*registeredStack.tiff']).name]), [4 2 6]);
            masks{jj} = tiffreadVolume([impath{ii}{jj} filesep 'maskStack.tiff']);
        end
    
        % Choose masks
        for jj = 1:numel(masks)
            [mask{jj}, note] = chooseMask(imgs{jj}, squeeze(masks{jj}));
            mask{jj} = sum(mask{jj}, 3); 
            if strcmp(note, 'F') || strcmp(note, 'X') 
                notes{end+1, 1} = impath{ii}{jj};
                notes{end, 2} = note;
            end
        end
    
        % Save masks
        for jj = 1:numel(mask)
            imwrite(mask{jj}, sprintf('%s%stotalMask.tiff', impath{ii}{jj}, filesep));
        end
    catch ME
        warning(getReport(ME))
        notes{end+1, 1} = impath{ii}{jj};
        notes{end, 3} = ME;
    end
end

% Save cell array with notes for later evalutaiton
save(sprintf('%s%smaskingNotes.mat', destination, filesep), 'notes');