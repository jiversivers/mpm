%% Find decay files
decay = FolderFinder('.sdt');
% The .sdt filenames default to the same name, but if they have been
% renamed, this line will hndle that.
fn = cellfun(@(x) dir([x filesep '*.sdt']).name, decay, 'UniformOutput', false);
decay_FilePath = cellfun(@(x, y) [x filesep y], decay, fn, 'UniformOutput', false);
treeRoot = findRoot(decay);
start = length(strsplit(treeRoot, filesep));
decayParts = pathpieces(decay);
%% Get calibration data with defaul parameters (will be asked for IRF file
% path)
calibration = IRFCalibration();
save('\\deckard\bmeg\Rajaram-Lab\Ivers,Jesse\Codes\MPM_Processing\FLIM Code\UprightIRF_1_calibration.mat', 'calibration');

%% Get G and S coords for all decays
G = cell(1, numel(decay_FilePath));
S = cell(size(G));

%% Run
for ii = 1:numel(decay_FilePath)
    fprintf('Processing %d/%d...\n', ii, numel(decay_FilePath))
    [G{ii}, S{ii}] = getPhasorCoords(decay_FilePath{ii}, calibration, 'harmonic', 1, 'triggerdelay', false, 'spatialbin', 2);
end
save([strrep(treeRoot, 'Raw', 'Processed') filesep 'GandSCoords_with_newIRF_10nsCollection.mat'], "G", "S", '-mat');
disp('Saved succesfully')

% Check errors
err = S(~cellfun(@isfloat, G));
% Looks like its all from the inverted and the files are too big to load
G = G(cellfun(@isfloat, G));
S = S(cellfun(@isfloat, S));
decay = decay(cellfun(@isfloat, G));
%% Process FOV multi-modal stuff
procRoot = '\\deckard\bmeg\Rajaram-Lab\Ivers,Jesse\UMSCC Radiation Study\MPM\Processed';
orrProc = FolderFinder('ORRmap.mat', procRoot);
pp = pathpieces(orrProc);
%%
errors = cell(0, 2);
findPartner = @(x, y) y(contains(y, x{end-2}) & contains(y, x{end-1}) & contains(y, x{end}));
for ii = 1:numel(orrProc)
    try
        partner = findPartner(pp(ii, :), decay_FilePath);
        if isempty(partner)
            continue
        end
    
        fprintf('\n%d/%d >>> %s\n\n', ii, numel(orrProc), orrProc{ii})
    
        % Calc phasor and load orr
        [g, s] = getPhasorCoords(partner{:}, calibration);
        load([orrProc{ii} filesep dir([orrProc{ii} filesep '*ORRmap.mat']).name], 'map');
        
        % Tau estimates mean and stack
        w = 2*pi*calibration.params.f;
        phi = atan(s./g);
        m = sqrt(g.^2 + s.^2);
        tau_phi = (1/w)*tan(phi);
        tau_m = (1/w)*sqrt((m.^-2) - 1);
    
        % Stack vectors
        ogs = cat(3, map, g, s);
        otptm = cat(3, map, g, s);
        
        % Save
        save([orrProc{ii} filesep 'orr_g_s.mat'], "ogs")
        save([orrProc{ii} filesep 'orr_tauphi_tauM.mat'], "otptm")
    catch ME
        errors{end+1, 1} = orrProc{ii};
        errors{end, 2} = getReport(ME);
        warning(getReport(ME))
    end

    % Try to add shg data
%     try
%         load([orrProc{ii} filesep dir([orrProc{ii} filesep '*_data_struct.mat']).name]);
%     catch
%         continue
%     end
%     
%     shg = Data([strcmp({Data.Wavelength}, '855')]).NormImg(:,:,4);
%     shgN = shg/max(shg, [], 'all');
%     col = shg.*(shgN > graythresh(shgN));
%     
end
    
%%
cat1 = unique(decayParts(:, start));
cat2 = unique(decayParts(:, start+1));
cat3 = unique(decayParts(:, start+2));

% Create and save cummulative phasors
w = 2*pi*calibration.params.f;
for ii = 1:numel(cat1)
    sub1 = [treeRoot filesep cat1{ii}];
    g = cat(3, G{strcmp(decayParts(:, start), cat1{ii})});
    s = cat(3, S{strcmp(decayParts(:, start), cat1{ii})});
    getAndSaveOutputs(g, s, w, ['\' strrep(sub1, 'Raw', 'Processed') filesep cat1{ii} '_freq10ns_'])
    for jj = 1:numel(cat2)
        sub2 = [sub1 filesep cat2{jj}];
        g = cat(3, G{strcmp(decayParts(:, start), cat1{ii}) & strcmp(decayParts(:, start+1), cat2{jj})});
        s = cat(3, S{strcmp(decayParts(:, start), cat1{ii}) & strcmp(decayParts(:, start+1), cat2{jj})});
        getAndSaveOutputs(g, s, w, ['\' strrep(sub2, 'Raw', 'Processed') filesep cat1{ii} '_' cat2{jj} '_freq10ns_'])
        for kk = 1:numel(cat3)
            sub3 = [sub2 filesep cat3{kk}];
            g = cat(3, G{strcmp(decayParts(:, start), cat1{ii}) & strcmp(decayParts(:, start+1), cat2{jj}) & strcmp(decayParts(:, start+2), cat3{kk})});
            s = cat(3, S{strcmp(decayParts(:, start), cat1{ii}) & strcmp(decayParts(:, start+1), cat2{jj}) & strcmp(decayParts(:, start+2), cat3{kk})});
            getAndSaveOutputs(g, s, w, ['\' strrep(sub3, 'Raw', 'Processed') filesep cat1{ii} '_' cat2{jj} '_' cat3{kk} '_freq10ns_'])
        end
    end
end

%%
IMProc = FolderFinder('orr_g_s.tiff');
pp = pathpieces(IMProc);
cat1 = unique(pp(:, 8));
cat2 = unique(pp(:, 9));
cat3 = unique(pp(:, 10));
%%
% Create and save cummulative phasors
w = 2*pi*calibration.params.f;
fig = figure;
fig.WindowState = 'fullscreen';
for ii = 1:numel(cat1)
    sub1 = [strrep(treeRoot, 'Raw', 'Processed') filesep cat1{ii}];
%     loadedOGS = cellfun(@(x) load([x filesep 'orr_g_s.tiff'], '-mat'), IMProc(strcmp(pp(:,8), cat1{ii})));
%     loadedOGS = squeeze(struct2cell(loadedOGS));
%     ogs = cat(3, loadedOGS{:});
%     o = ogs(:,:,1:3:end);
%     g = ogs(:,:,2:3:end);
%     s = ogs(:,:,3:3:end);
%     phi = atan(s./g);
%     m = sqrt(g.^2 + s.^2);
%     tau_phi = (1/w)*tan(phi);
%     tau_m = (1/w)*sqrt((m.^-2) - 1);
    
    for jj = 1:numel(cat2)
        sub2 = [sub1 filesep cat2{jj}];
        for kk = 1:numel(cat3)
            sub3 = [sub2 filesep cat3{kk}];
            loadedOGS = cellfun(@(x) load([x filesep 'orr_g_s.tiff'], '-mat'), IMProc(strcmp(pp(:,8), cat1{ii}) & strcmp(pp(:,9), cat2{jj}) & strcmp(pp(:,10), cat3{kk})));
            loadedOGS = squeeze(struct2cell(loadedOGS));
            ogs = cat(3, loadedOGS{:});
            o = ogs(:,:,1:3:end);
            g = ogs(:,:,2:3:end);
            s = ogs(:,:,3:3:end);
            phi = atan(s./g);
            m = sqrt(g.^2 + s.^2);
            tau_phi = (1/w)*tan(phi);
            tau_m = (1/w)*sqrt((m.^-2) - 1); 
            scatter3(o(:), g(:), s(:), '.k');
            ax = gca;
            ax.XLabel.String = 'ORR';
            ax.YLabel.String = 'G';
            ax.ZLabel.String = 'S';
            ax.XLabel.FontSize = 24;
            ax.YLabel.FontSize = 24;
            ax.ZLabel.FontSize = 24;
            exportgraphics(gcf, [sub3 filesep 'ogs_3dScatter.png']);
            scatter(o(:), g(:), '.k');
            ax.XLabel.String = 'ORR';
            ax.YLabel.String = 'G';
            ax.XLabel.FontSize = 24;
            ax.YLabel.FontSize = 24;
            exportgraphics(gcf, [sub3 filesep 'og_2dScatter.png']);
            scatter(o(:), s(:), '.k');
            ax.XLabel.String = 'ORR';
            ax.YLabel.String = 'S';
            ax.XLabel.FontSize = 24;
            ax.YLabel.FontSize = 24;
            exportgraphics(gcf, [sub3 filesep 'os_2dScatter.png']);
        end
    end
end

%%
function getAndSaveOutputs(g, s, w, savename)
    phi = atan(s./g);
    m = sqrt(g.^2 + s.^2);
    tau_phi = mean((1/w)*tan(phi), 'all', 'omitnan');
    tau_m = mean((1/w)*sqrt((m.^-2) - 1), 'all', 'omitnan');
    writecell({tau_phi, tau_m}, [savename '_lifetimemeans.xlsx']);
    plotPhasor(g(:), s(:));
    phasorView
    ax = gca;
    ax.XLabel.String = 'G'; ax.YLabel.String= 'S';
    cbar = colorbar;
    cbar.Label.String = 'Normalized Frequency';
    fig = gcf;
    fig.WindowState = 'fullscreen';
    exportgraphics(fig, [savename '_phasor.png'])
end

