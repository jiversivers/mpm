%% Bulk file handling
% Get all files and determine types
startDir = uigetdir(pwd, 'Select data root...');
pwrDir = uigetdir(startDir, 'Select power directory...');
saveDir = uigetdir(startDir, 'Select save root...');
if ~exist([saveDir filesep 'Processed'], 'dir')
    mkdir([saveDir filesep 'Processed'])
end
[irfFile, irfPath] = uigetfile({'*.sdt', 'Raw Decays'; '*tif; *.tiff', 'Image Stacks'; '*.mat', 'Calibration Files'}, 'Select IRF Calibration file...', startDir);
[~, ~, ext] = fileparts(irfFile);

%% Parse and process ui-inputs
if any(strcmp(ext, {'.sdt', '.tif', 'tiff'}))
    params = IRFCalibration([irfPath irfFile]);
else
    load(irfFile, 'params')
end
rawDir = FolderFinder('References', startDir);
flim = rawDir(contains(rawDir, 'flim', 'IgnoreCase', true));
rdox = rawDir(contains(rawDir, 'redox','IgnoreCase', true));

fits = FolderFinder('all.asc', saveDir);

% Back img dirs out to find common FOV image and group dirs
allParts = pathpieces(horzcat(flim, rdox));
fov = unique(namemaker(allParts(:,1:end-1), filesep));
group = unique(namemaker(allParts(:,1:end-5), filesep));

% Get root for processed file saving
sd = pathpieces(startDir);

% Prep output table
catHdrs = 'Cell Line\\Treatment\\Hours Post Treatment\\Animal ID\\Sample\\FOV';
varHdrs = '\\%1$s FAD\\%1$s NADH\\%1$s ORR from Map\\%1$s SHG\\%1$s G\\%1$s S\\%1$s Photons\\%1$s TauPhi\\%1$s TauM\\%1$s Mean Tau\\%1$s A2Percent\\%1$s Chi2\\%1$s Tm Good-Fits\\%1$s A2Percent Good-Fits';
hdrs = strsplit(horzcat(catHdrs, sprintf(varHdrs, 'Avg'), sprintf(varHdrs, 'Var')), '\');
outTab = cell(1+numel(fov), numel(hdrs));
outTab(1,:) = hdrs;

%% Process all FOV
warning('off', 'MATLAB:imagesci:Tiff:missingExtraSamples');
p = 512^2;
cmap = colormap('jet').*colormap('gray');
r = 1; % Start at 1 to handle header row
% Group-level cummulative outputs
for ii = 1:numel(group)
    gr = group(ii);
    FOV = fov(contains(fov, gr));
    N = numel(FOV);

    % These variable will accumulate the group-wide FOV variables of the
    % same name (lowercase)
    ORR = NaN(N*p, 1);
    ORRMask = NaN(N*p, 1);
    SHG = NaN(N*p, 1);
    G = NaN(N*p, 1);
    S = NaN(N*p, 1);
    Photons = NaN(N*p, 1);
    fitG = NaN(N*p, 1);
    fitS = NaN(N*p, 1);
    fitPhotons = NaN(N*p, 1);
    Tm = NaN(N*p, 1);
    A2Percent = NaN(N*p, 1);
    X2 = NaN(N*p, 1);

    % Process all individual FOV outputs and accumulate group's vars for
    % group porcessing of all FOVs w/in group
    for jj = 1:numel(FOV)
        fprintf('Group %d, FOV %d of %d.', ii, jj, numel(FOV))
        f = FOV(jj);
        r = r+1;

        %{
        Try/catch blocks will handle cases where flim or redox data wasn't
        obtained.
        %}

        % Redox processing
        try
            img = getAFSingalsFromDirs(rdox{contains(rdox, f)}, pwrDir);
            if size(img.NADH) ~= [512 512]
                warning('Resizing %s\n', f{:}')
                img = structfun(@(x) imresize(x, [512 512 []]), img, 'UniformOutput', false);
            end
        catch ME
            warning('%s\n occured when processing redox for %s.\n', getReport(ME), f{:});
            img.FAD = NaN(512, 512);
            img.NADH = NaN(512, 512);
            img.SHG = NaN(512, 512);
            img.all = NaN(512, 512, 8);
            img.n = 0;
        end
        orrMask = (img.FAD > 3e-4 & img.NADH > 3e-4) & ~any(isnan(img.all), 3);
        img.n = sum(orrMask, 'all');
        orr = img.FAD./(img.FAD+img.NADH);
        orr(~orrMask) = nan;
        int = (img.FAD+img.NADH)/2;
        int(~orrMask) = nan;

        % Phasor Processing
        try
            sdt = [flim{contains(flim, f)} filesep dir([flim{contains(flim, f)} filesep '*.sdt']).name];
            [g, s, ph] = getPhasorCoords(sdt, params);
            N = numel(ph);
        catch ME
            warning('%s\n occured when processing FLIM for %s.\n', getReport(ME), f{:});
            g = NaN(512, 512);
            s = NaN(512, 512);
            ph = NaN(512, 512);
            N = 0;
        end
        phasOut = getPhasorOutputs(g, s, params.w);

        % Fits Processing
        try
            fitData = loadFitData(fits{contains(fits, strrep(f, 'Raw', 'FLIM Fit Outs'))});
            fitData.n = numel(fitData.photons);
        catch ME
            warning('%s\n occured when loading fits for %s.\n', getReport(ME), f{:});
            fitData.G = NaN(512, 512);
            fitData.S = NaN(512, 512);
            fitData.tm = NaN(512, 512);
            fitData.a2percent = NaN(512, 512);
            fitData.photons = NaN(512, 512);
            fitData.chi = 1000000*ones(512, 512);
            fitData.n = 0;
        end
        goodFit = isWithin(fitData.chi, [0 1.4], '=');

        % Append FOV stats to table
        avgOuts = cellfun(@(x) mean(x, 'all', 'omitnan'), {img.FAD, img.NADH, orr, img.SHG, g, s, ph, phasOut.Tau.Phi, phasOut.Tau.M, fitData.tm, fitData.a2percent, fitData.chi, fitData.tm(goodFit), fitData.a2percent(goodFit)});
        varOuts = cellfun(@(x) var(x, 0, 'all', 'omitnan'), {img.FAD, img.NADH, orr, img.SHG, g, s, ph, phasOut.Tau.Phi, phasOut.Tau.M, fitData.tm, fitData.a2percent, fitData.chi, fitData.tm(goodFit), fitData.a2percent(goodFit)});
        nOuts = {img.n, N, fitData.n};
        fp = pathpieces(f);
        outTab(r, 1:6) = fp(1, [4:end-2, end]);
        outTab(r, 7:end) = num2cell([avgOuts, varOuts]);

        % Update cummulative vars
        ORR(p*(jj-1)+1:p*jj, :) = orr(:);
        ORRMask(p*(jj-1)+1:p*jj, :) = orrMask(:);
        SHG(p*(jj-1)+1:p*jj, :) = img.SHG(:);
        G(p*(jj-1)+1:p*jj, :) = g(:);
        S(p*(jj-1)+1:p*jj, :) = s(:);
        Photons(p*(jj-1)+1:p*jj, :) = ph(:);
        fitG(p*(jj-1)+1:p*jj, :) = fitData.G(:);
        fitS(p*(jj-1)+1:p*jj, :) = fitData.S(:);
        fitPhotons(p*(jj-1)+1:p*jj, :) = fitData.photons(:);
        Tm(p*(jj-1)+1:p*jj, :) = fitData.tm(:);
        A2Percent(p*(jj-1)+1:p*jj, :) = fitData.a2percent(:);
        X2(p*(jj-1)+1:p*jj, :) = fitData.chi(:);
        
        % Write FOV outputs
        fovSave = namemaker(horzcat({saveDir, 'Processed'}, fp(1,numel(sd)+1:end)), filesep);
        if ~exist(fovSave, 'dir')
            mkdir(fovSave)
        end

        try 
            TiffSingleWrite(cat(3, img.all(:,:,5), img.FAD, img.NADH, orr, img.SHG, g, s, phasOut.Tau.Phi, phasOut.Tau.M, fitData.tm, fitData.a2percent, fitData.chi), ...
                            [fovSave filesep 'outputMapStack.tiff'], 1);
            orr_fig = createColorORRMap(orr, int, 'redoxMax', 0.8);
            imwrite(orr_fig, [fovSave filesep 'coloredORR.tiff'])

        catch ME
            warning(getReport(ME))
        end
    end
    
    % Save group wide vars for future use
    gp = pathpieces(gr);
    grpSave = namemaker(horzcat({saveDir, 'Processed'}, gp(1,numel(sd)+1:end)), filesep);
        if ~exist(grpSave, 'dir')
            mkdir(grpSave)
        end
    save([grpSave filesep 'ORR_Mask_SHG_G_S_Ph_Tm_A2Percent_X2.mat'], "ORR", "ORRMask", "SHG", "G", "S", "Photons", "Tm", "A2Percent", "X2", '-mat');

    % Reshape vars back into image sequences (preparing for 2D medfilts)
    G = reshape(G, 512, []);
    S = reshape(S, 512, []);
    Photons = reshape(Photons, 512, []);
    
    fitG = reshape(fitG, 512, []);
    fitS = reshape(fitS, 512, []);
    fitPhotons = reshape(fitPhotons, 512, []);
    
    ORR = reshape(ORR, 512, []);

    Tm = reshape(Tm, 512, []);
    A2Percent = reshape(A2Percent, 512, []);
    X2 = reshape(X2, 512, []);

    % Create cummulative group plots and save
    % Set thresholds/limits
    chiLim = [0 1.4];
    wInChiLims = @(x) x>=chiLim(1) & x<=chiLim(2);
    Tm(~wInChiLims(X2)) = NaN;
    A2Percent(~wInChiLims(X2)) = NaN;
    ptTh = 15;

    % First, my phasors:
    % Raw phasor plot
    fig = figure;
    ax = axes(fig);
    phasor = plotPhasor(ax, G, S, 'photons', Photons, 'intensitythreshold', 0, 'nbins', [128 128]);
    phasorView(phasor, cmap, [0.5 1 2 4 8], params)
    exportgraphics(phasor.Parent, [grpSave filesep 'myRawPhasor.png'], 'Resolution', 300)
    savefig(gcf, [grpSave filesep 'myRawPhasor.fig'])

    % Phasor plot with single pass [3x3] medfilt
    fig = figure;
    ax = axes(fig);
    phasor = plotPhasor(ax, G, S, 'photons', Photons, 'intensitythreshold', 0, 'nbins', [128 128], 'medfiltsize', [3 3], 'medfiltcount', 1);
    phasorView(phasor, cmap, [0.5 1 2 4 8], params)
    exportgraphics(phasor.Parent, [grpSave filesep 'myMedfilt3x3Phasor.png'], 'Resolution', 300)
    savefig(gcf, [grpSave filesep 'myMedfilt3x3Phasor.fig'])

    % Phasor plot with single pass [3x3] medfilt and intensity threshold
    % of phTh
    phTh = 15;
    fig = figure;
    ax = axes(fig);
    [phasor, N, ctrs] = plotPhasor(ax, G, S, 'photons', Photons, 'intensitythreshold', phTh, 'nbins', [128 128], 'medfiltsize', [3 3], 'medfiltcount', 1);
    phasorView(phasor, cmap, [0.5 1 2 4 8], params)
    exportgraphics(phasor.Parent, [grpSave filesep 'myMedfilt3x3PhasorIntThresh' num2str(phTh) '.png'], 'Resolution', 300)
    savefig(gcf, [grpSave filesep 'myMedfilt3x3PhasorIntThresh' num2str(phTh) '.fig'])

    % Phasor (x,y,z) v ORR (color) w/ medfilt and int thresh
    fig = figure;
    ax = axes(fig);
    
    % Medfilt vars for 3D Histogrm
    Var1 = medfilt2(G, [3 3]);
    Var2 =  medfilt2(S, [3 3]);
    colorBasis = ORR;
    
    % Threshold by photon counts
    Var1(Photons<phTh) = NaN;
    Var2(Photons<phTh) = NaN;
    colorBasis(Photons<phTh) = NaN;
    
    % Create and plot 3d ColoredHist of G v S v ORR
    C = HistogramColorCalculator(Var1, Var2, colorBasis, ctrs);
    surf(ax, ctrs{:}, N, C)
    shading interp
    view(3)
    ax.XLabel.String = 'G';
    ax.YLabel.String = 'S';
    ax.ZLabel.String = '% of Total Pixels';
    ax.XLim = [0 1];
    ax.YLim = [0 0.75];
    ax.ZLim = [0 0.75];
    ax.XTick = 0:0.5:1;
    ax.YTick = [0 0.5];
    ax.ZTick = 0:.25:.75;
    pbaspect([4 3 1])
    viscircles(ax, [0.5, 0], 0.5);
    colormap(cmap);
    clim([0 0.8])
    cbar = colorbar;
    cbar.Label.String = 'ORR';
    exportgraphics(ax, [grpSave filesep 'my3d_Redox_v_3x3MedFilt_PhotonThresh' num2str(phTh) '_Lifetime_Histogram.png'], 'Resolution', 300)
    savefig(gcf, [grpSave filesep 'my3d_Redox_v_3x3MedFilt_PhotonThresh' num2str(phTh) '_Lifetime_Histogram.fig']);

    % Now, SPCImage phasors:
    % Raw phasor plot
%     fig = figure;
%     ax = axes(fig);
%     phasor = plotPhasor(ax, fitG, fitS, 'photons', fitPhotons, 'intensitythreshold', 0, 'nbins', [128 128]);
%     phasorView(phasor, cmap, [0.5 1 2 4 8], params)
%     exportgraphics(phasor.Parent, [grpSave filesep 'spcRawPhasor.png'], 'Resolution', 300)
%     savefig(gcf, [grpSave filesep 'spcRawPhasor.fig'])
% 
%     % Phasor plot with single pass [3x3] medfilt
%     fig = figure;
%     ax = axes(fig);
%     phasor = plotPhasor(ax, fitG, fitS, 'photons', fitPhotons, 'intensitythreshold', 0, 'nbins', [128 128], 'medfiltsize', [3 3], 'medfiltcount', 1);
%     phasorView(phasor, cmap, [0.5 1 2 4 8], params)
%     exportgraphics(phasor.Parent, [grpSave filesep 'spcMedfilt3x3Phasor.png'], 'Resolution', 300)
%     savefig(gcf, [grpSave filesep 'spcMedfilt3x3Phasor.fig'])
% 
%     % Phasor plot with single pass [3x3] medfilt and intensity threshold
%     % of phTh
%     fig = figure;
%     ax = axes(fig);
%     [phasor, N, ctrs] = plotPhasor(ax, fitG, fitS, 'photons', fitPhotons, 'intensitythreshold', phTh, 'nbins', [128 128], 'medfiltsize', [3 3], 'medfiltcount', 1);
%     phasorView(phasor, cmap, [0.5 1 2 4 8], params)
%     exportgraphics(phasor.Parent, [grpSave filesep 'spcMedfilt3x3PhasorIntThresh' num2str(phTh) '.png'], 'Resolution', 300)
%     savefig(gcf, [grpSave filesep 'scpMedfilt3x3PhasorIntThresh' num2str(phTh) '.fig'])
%     
%     % Phasor (x,y,z) v ORR (color) w/ medfilt and int thresh
%     fig = figure;
%     ax = axes(fig);
%     
%     % Medfilt vars for 3D Histogrm
%     Var1 = medfilt2(fitG, [3 3]);
%     Var2 =  medfilt2(fitS, [3 3]);
%     colorBasis = ORR;
%     
%     % Threshold by photon counts
%     Var1(fitPhotons<phTh) = NaN;
%     Var2(fitPhotons<phTh) = NaN;
%     colorBasis(fitPhotons<phTh) = NaN;
%     
%     % Create and plot 3d ColoredHist of G v S v ORR
%     C = HistogramColorCalculator(Var1, Var2, colorBasis, ctrs);
%     surf(ax, ctrs{:}, N, C)
%     shading interp
%     view(3)
%     ax.XLabel.String = 'G';
%     ax.YLabel.String = 'S';
%     ax.ZLabel.String = '% of Total Pixels';
%     ax.XLim = [0 1];
%     ax.YLim = [0 0.75];
%     ax.ZLim = [0 1.1];
%     ax.XTick = 0:0.5:1;
%     ax.YTick = [0 0.5];
%     ax.ZTick = [0 .5 1];
%     pbaspect([4 3 1])
%     viscircles(ax, [0.5, 0], 0.5);
%     colormap(cmap);
%     clim([0 0.8])
%     cbar = colorbar;
%     cbar.Label.String = 'ORR';
%     exportgraphics(ax, [grpSave filesep 'spc3d_Redox_v_3x3MedFilt_PhotonThresh' num2str(phTh) '_Lifetime_Histogram.png'], 'Resolution', 300)
%     savefig(gcf, [grpSave filesep 'spc3d_Redox_v_3x3MedFilt_PhotonThresh' num2str(phTh) '_Lifetime_Histogram.fig']);

    % Fit-based Histograms
    fig = figure;
    ax = axes(fig);
    
    % Vars for 3D Histogrm
    Var1 = Tm;
    Var2 =  A2Percent;
    colorBasis = ORR;
    
    % Threshold by fit quality
    colorBasis(~wInChiLims(X2)) = NaN;

    % 2D hist of Tm v A2Per
    ctrs{1} = linspace(0, 3500, 128);
    ctrs{2} = linspace(0, 100, 128);
    N = hist3([Var1(:), Var2(:)], 'Ctrs', ctrs);
    N = 100*(N'/sum(N, 'all'));
    surf(ax, ctrs{:}, N);
    view(2)
    shading interp
    ax.XLabel.String = '\tau_M (ps)';
    ax.YLabel.String = '\alpha_2%';
    colormap(cmap)
    clim([0 .85])
    cbar = colorbar;
    cbar.Label.String = '% of Total Pixels';
    exportgraphics(ax, [grpSave filesep 'spc2d_TmvA2Percent_Histogram.png'], 'Resolution', 300)
    savefig(gcf, [grpSave filesep 'spc2d_TmvA2Percent_Histogram.fig']);

    % Fit-based 3d Histogram
    fig = figure;
    ax = axes(fig);
    
    % Create and plot 3d ColoredHist of Tm v A2 v ORR
    C = HistogramColorCalculator(Var1, Var2, colorBasis, ctrs);
    surf(ax, ctrs{:}, N, C)
    shading interp
    view(3)
    ax.XLabel.String = '\tau_M (ps)';
    ax.YLabel.String = '\alpha_2%';
    ax.ZLabel.String = '% of Total Pixels';
    ax.XLim = [0 3500];
    ax.YLim = [0 100];
    ax.ZLim = [0 75];
    ax.XTick = 0:500:3500;
    ax.YTick = 0:25:100;
    ax.ZTick = 0:25:75;
    colormap(cmap);
    clim([0 0.8])
    cbar = colorbar;
    cbar.Label.String = 'ORR';
    exportgraphics(ax, [grpSave filesep 'spc3d_Redox_v_Fit__Lifetime_Histogram.png'], 'Resolution', 300)
    savefig(gcf, [grpSave filesep 'spc3d_Redox_v_Fit__Lifetime_Histogram.fig']);

    % 1D Hist of ORR
    fig = figure;
    ax = axes(fig);
    [N, edges] = histcounts(ORR, 25, 'BinLimits', [0 1], 'Normalization', 'probability');
    plot(ax, (edges(3:end-1)+edges(2:end-2))/2, 100*N(2:end-1), 'LineWidth', 2.5)
    ax.XLabel.String = 'ORR';
    ax.YLabel.String = '% of Total Pixels in Group';
    ax.XLim = [0 1];
    ax.YLim = [0 10];
    ax.YScale = 'linear';
    exportgraphics(ax, [grpSave filesep 'orr_histogram.png'], 'Resolution', 300)
    savefig(gcf, [grpSave filesep 'orr_histogram.fig']);

    % 1D Hist of Tm
    fig = figure;
    ax = axes(fig);
    [N, edges] = histcounts(Tm, 25, 'BinLimits', [0 3500], 'Normalization', 'probability');
    plot(ax, (edges(3:end-1)+edges(2:end-2))/2, 100*N(2:end-1), 'LineWidth', 2.5)
    ax.XLabel.String = '\tau_M (ps)';
    ax.YLabel.String = '% of Total Pixels in Group';
    ax.XLim = [0 3500];
    ax.YLim = [0 25];
    ax.YScale = 'linear';
    exportgraphics(ax, [grpSave filesep 'tm_histogram.png'], 'Resolution', 300)
    savefig(gcf, [grpSave filesep 'tm_histogram.fig']);

    % 1D Hist of A2Percent
    fig = figure;
    ax = axes(fig);
    [N, edges] = histcounts(A2Percent, 25, 'BinLimits', [0 100], 'Normalization', 'probability');
    plot(ax, (edges(3:end-1)+edges(2:end-2))/2, 100*N(2:end-1), 'LineWidth', 2.5)
    ax.XLabel.String = '\alpha_2%';
    ax.YLabel.String = '% of Total Pixels in Group';
    ax.XLim = [0 100];
    ax.YLim = [0 50];
    ax.YScale = 'linear';
    exportgraphics(ax, [grpSave filesep 'a2percent_histogram.png'], 'Resolution', 300)
    savefig(gcf, [grpSave filesep 'a2percent_histogram.fig']);

    close all
end
writecell(outTab, 'fov_outputs.xlsx');
warning('on', 'MATLAB:imagesci:Tiff:missingExtraSamples');