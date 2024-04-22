% Find all group terminals
terms = FolderFinder('sam');
nP = pathpieces(terms);

%%
% Go there and make the cummPlot
fN = cell(0);
for ii = 1:length(terms)
    [fig, fN{ii}] = cummHistPlot(terms{ii}, '-orr');
end

%% Make Plot of all
figure;
h2 = tiledlayout(4,3);
for ii = 1:length(terms)
    h1 = openfig(fN{ii}); pause(1)
    ax = gca;pause(1)
    ax.Parent = h2;
    ax.Layout.Tile = ii;
    
    title([nP{ii, end-1}, ' ', nP{ii, end}])
    close(h1)
end

%%
function [fig, filename] = cummHistPlot(startDir, typeFlag)
    %% Gather Data
    cd(startDir)

    if strcmp(typeFlag, '-orr')
        keyword = '_ORRmap.mat';
    
        % Find all ORR maps in the group
        fov = FolderFinder(keyword);
        
        % Combine all ORR vals into one vector
        groupVals = [];
        groupWeights = [];
        jj = 1;
        for ii = 1:length(fov)
            cd(fov{ii})
        
            % Load variables
            load(dir('*_ORRmap.mat').name);
            load(dir('*_data_struct.mat').name);
            
            % Generate intensity map (avg of NADH and FAD intensities)
            for ii = 1:length(Data)
                if Data(ii).Wavelength == '755'
                    nadh = Data(ii).UseImg(:,:,3);
                elseif Data(ii).Wavelength == '855'
                    fad = Data(ii).UseImg(:,:,2);
                end
            end
            
            int = (nadh+fad)./2;
            
        
            % Ignore background and NaN
            vals =  map(map~=0 & ~isnan(map));
            wats = int(map~=0 & ~isnan(map));
            sz = numel(vals);
            groupVals(jj:sz+jj-1) = vals;
            groupWeights(jj:sz+jj-1) = wats;
            jj = jj + sz;
        end

    elseif strcmp(typeFlag, '-taum')
        keyword = '_t1.tiff';
        
        % Find all maps in the group
        fov = FolderFinder(keyword);
        
        % Combine all vals into one vector
        groupVals = [];
        groupWeights = [];
        jj = 1;
        for ii = 1:length(fov)
            cd(fov{ii})
        
            % Load images and calculate TauM
            tau1 = double(imread(dir('*_t1.tiff').name));
            tau2 = double(imread(dir('*_t2.tiff').name));
            alp1 = double(imread(dir('*_a1[%].tiff').name));
            alp2 = double(imread(dir('*_a2[%].tiff').name));
            map = (alp1/100).*tau1 + (alp2/100).*tau2;

            int = double(imread(dir('*_intensity_image.tif').name));

            % Ignore background and NaN
            vals =  map(map>=0 & ~isnan(map));
            wats = int(map>=0 & ~isnan(map));
            sz = numel(vals);
            groupVals(jj:sz+jj-1) = vals;
            groupWeights(jj:sz+jj-1) = wats;
            jj = jj + sz;
        end

    elseif strcmp(typeFlag, '-a1p')
        keyword = '_a1[%].tiff';
        
        % Find all maps in the group
        fov = FolderFinder(keyword);
        
        % Combine all vals into one vector
        groupVals = [];
        groupWeights = [];
        jj = 1;
        for ii = 1:length(fov)
            cd(fov{ii})
        
            % Load images
            map = double(imread(dir('*_a1[%].tiff').name))/100;
            int = double(imread(dir('*_intensity_image.tif').name));

            % Ignore background and NaN
            vals =  map(map>=0 & ~isnan(map));
            wats = int(map>=0 & ~isnan(map));
            sz = numel(vals);
            groupVals(jj:sz+jj-1) = vals;
            groupWeights(jj:sz+jj-1) = wats;
            jj = jj + sz;
        end
    else
        error('Unrecognized type flag. Use -orr, -taum, or -a1p.')
    end
    
    %% Make fig
    
    % Prepare histogram counts
    [N, centers] = weightedHistogram(groupVals, groupWeights, 100);
    N = N(2:end-1);
    N = N/max(N);
    centers = centers(2:end-1);
    
    % Plot histogram data
    br = bar(centers, N,'hist');
    br.FaceColor = [0.75 0.75 0.75];
    br.EdgeColor = [0.5 0.5 0.5];
    br.FaceAlpha = 0.5;
    br.EdgeAlpha = 0.5;
    
    % Scale figure
    fig = gcf;
    ax = gca;
    if strcmp(typeFlag, '-orr')
        ax.YLim = [0 1];
        ax.XLim = [0 1];
    elseif strcmp(typeFlag, '-taum')
        ax.YLim = [0 1];
        ax.XLim = [0 4000];
    elseif strcmp(typeFlag, '-a1p')
        ax.YLim = [0 1];
        ax.XLim = [0 1000];
    end

    
    % Save individual figure
    cd(startDir)
    nP = pathpieces(pwd);
    filename = [startDir, filesep, nP{end-1},'_',nP{end}, typeFlag];
    disp(['Saving figure for ', filename])
    savefig(fig, filename)
end
