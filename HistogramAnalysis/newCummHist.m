root = pwd;
histData = struct();

orrEdges = [0:0.02:1];
t2Edges = [0:1000:100000];
A1A2Edges = [0:0.05:5];

ctrsG = 0+0.5*0.01:0.01:1-0.5*0.01;
ctrsS = 0+0.5*0.005:0.005:0.5-0.5*0.005;
PhasCtrs = {ctrsG, ctrsS};

t = table();


for Cl = {'UMSCC22B','UMSCC47'}
    for Tx = {'NT', 'XT'}
        for Ht = {'Baseline', '1hpt', '24hpt', '48hpt'}
            currentDir = [root, filesep, Cl{:}, filesep, Tx{:}, filesep, Ht{:}];
            try
                cd(currentDir)
            catch
                continue
            end
            terms = FolderFinder('_data_struct.mat');
            orrHist = zeros(1,100);
            orrVals = [];
            intVals = [];
            GVals = [];
            SVals = [];
            t2Hist = zeros(1,100);
            A1A2Hist = zeros(1,100);
            PhasHist = zeros(100, 100);
            for ii = 1:numel(terms)
                try
                % Load images
                load([terms{ii}, filesep, dir([terms{ii}, filesep, 'fov*_data_struct.mat']).name])
                phot = double(imread([terms{ii}, filesep, dir([terms{ii}, filesep, '*_photons.tiff']).name]));
                t1 = imread([terms{ii}, filesep, dir([terms{ii}, filesep, '*_t1.tiff']).name]);
                t2 = imread([terms{ii}, filesep, dir([terms{ii}, filesep, '*_t2.tiff']).name]);
                a1 = imread([terms{ii}, filesep, dir([terms{ii}, filesep, '*_a1.tiff']).name]);
                a2 = imread([terms{ii}, filesep, dir([terms{ii}, filesep, '*_a2.tiff']).name]);
                
                phG = double(imread([terms{ii}, filesep, dir([terms{ii}, filesep, '*_phasor_G.tiff']).name]));
                phS = double(imread([terms{ii}, filesep, dir([terms{ii}, filesep, '*_phasor_S.tiff']).name]));
                
                % Fix my normalization mistakes
                % For NADH
                is755 = strcmp({Data.Wavelength}, '755');
                tempNadh = uint16(Data(is755).RawImg(:,:,3) + Data(is755).Modes(3));
                gan = Data(is755).Gain(3);
                pok = Data(is755).Pockels;
                rawPmp = keys(Data(is755).Power.Laser);
                rawLas = values(Data(is755).Power.Laser);
                rawObj = values(Data(is755).Power.Objective);
                mis = cellfun(@ismissing, rawPmp) + cellfun(@ismissing, rawLas) + cellfun(@ismissing, rawObj);
                pmp = cell2mat(rawPmp(~mis));
                las = cell2mat(rawLas(~mis));
                obj = cell2mat(rawObj(~mis));
                dte = Data(is755).Date;
                        % Get laser to add flag for normalization
                if Data(is755).Laser == "MaiTai"
                    lasFlag = '-m';
                elseif Data(is755).Laser == "InsightX3"
                    lasFlag = '-i';
                else
                    warning('Unrecognized laser.')
                    lasFlag = '-u';
                end
                nadh = fluorescein_normalization_master(tempNadh, gan, pok, '-las', pmp, las, obj, dte, lasFlag);
                
                % for FAD
                is855 = strcmp({Data.Wavelength}, '855');
                tempFad = uint16(Data(is855).RawImg(:,:,2) + Data(is855).Modes(2));
                gan = Data(is855).Gain(2);
                pok = Data(is855).Pockels;
                rawPmp = keys(Data(is855).Power.Laser);
                rawLas = values(Data(is855).Power.Laser);
                rawObj = values(Data(is855).Power.Objective);
                mis = cellfun(@ismissing, rawPmp) + cellfun(@ismissing, rawLas) + cellfun(@ismissing, rawObj);
                pmp = cell2mat(rawPmp(~mis));
                las = cell2mat(rawLas(~mis));
                obj = cell2mat(rawObj(~mis));
                fad = fluorescein_normalization_master(tempFad, gan, pok, '-las', pmp, las, obj, dte, lasFlag);

                %
                int = (double(nadh)+double(fad))/2;
                map = double(fad)./(double(nadh)+double(fad));
                mask = ~isnan(map);

                % Bulk ORR
                bNadh = sum(nadh, 'all');
                bFad = sum(fad, 'all');
                bORR = bFad/(bFad+bNadh);

                % Pixel ORR
                pORR = mean(map(mask));
                nP = pathpieces(terms{ii});
                save([terms{ii}, filesep, nP{end}, '_ORRmap.mat'], "map");

                % Lifetime stats
                tm = a1.*t1 + a2.*t2;
                mtm = mean(tm, 'all');
                A1A2 = a1./a2;
                mA1A2 = mean(A1A2, 'all', 'omitnan');

                % Table
                t(end+1, 1:numel(nP)) = nP;
                t(end, {'Bulk ORR', 'Pixel-wise ORR', 'Mean Lifetime', 'A1/A2 Ratio'}) = {bORR, pORR, mtm, mA1A2};


                % Nanmask and hist each
                orrVals = cat(1, orrVals, map(mask));
                intVals = cat(1, intVals, int(mask));
                N = histcounts(map(mask), 100);
                orrHist = orrHist+N;

%                 t2Vals = t2(mask & ~isnan(t2) & map~=0);
%                 N = histcounts(t2Vals, t2Edges);
%                 t2Hist = t2Hist+N;
% 
%                 A1A2Vals = A1A2(mask & ~isnan(A1A2) & map~=0);
%                 N = histcounts(A1A2Vals, A1A2Edges);
%                 A1A2Hist = A1A2Hist+N;

                mask = phot >= 120;
                badG = phG==0 | phG>1;
                badS = phS==0 | phS>1;
                GVals = cat(1, GVals, phG(mask & ~badG & ~badS));
                SVals = cat(1, SVals, phS(mask & ~badG & ~badS));
                N = hist3([phG(mask & ~badG & ~badS), phS(mask & ~badG & ~badS)], 'Ctrs', PhasCtrs);
                PhasHist = PhasHist + N;
                catch ME
                    warning(['Error in ', terms{ii}])
                end
            end
            % Normalize to pixel fraction
            orrHist = orrHist/sum(orrHist);
            t2Hist = t2Hist/sum(t2Hist);
            A1A2Hist = A1A2Hist/sum(A1A2Hist);
            PhasHist = PhasHist/sum(PhasHist, 'all');

            histData.(Cl{:}).(Tx{:}).(matlab.lang.makeValidName(Ht{:})).ORR = orrHist;
            histData.(Cl{:}).(Tx{:}).(matlab.lang.makeValidName(Ht{:})).ORRRaw = orrVals;
            histData.(Cl{:}).(Tx{:}).(matlab.lang.makeValidName(Ht{:})).Weights = intVals;
%             histData.(Cl{:}).(Tx{:}).(matlab.lang.makeValidName(Ht{:})).Tau2 = t2Hist;
%             histData.(Cl{:}).(Tx{:}).(matlab.lang.makeValidName(Ht{:})).A1A2 = A1A2Hist;
            histData.(Cl{:}).(Tx{:}).(matlab.lang.makeValidName(Ht{:})).Phas = PhasHist;
            histData.(Cl{:}).(Tx{:}).(matlab.lang.makeValidName(Ht{:})).GRaw = GVals;
            histData.(Cl{:}).(Tx{:}).(matlab.lang.makeValidName(Ht{:})).SRaw = SVals;

        end
    end
end

%%
cntrs = (orrEdges(2:end) + orrEdges(1:end-1))/2;
for Cl= {'UMSCC22B','UMSCC47'}
    for Tx = {'NT', 'XT'}
            fig = figure; hold on
            fig.Position = [1 41 1600 783];
            title([Cl{:}, ' ', Tx{:}], 'FontSize',30)
            xlabel('ORR', 'FontSize',20)
            ylabel('Pixel Fraction', 'FontSize',20)
        for Ht = {'Baseline',  '24hpt', '48hpt'}
            try
                plot(cntrs, histData.(Cl{:}).(Tx{:}).(matlab.lang.makeValidName(Ht{:})).ORR, 'LineWidth', 5)           
            catch
                plot(cntrs, histData.(Cl{:}).NT.(matlab.lang.makeValidName(Ht{:})).ORR, 'LineWidth', 5)
                continue
            end
        end
        ylim([0, 0.05])
        legend({'Baseline',  '24hpt', '48hpt'}, 'FontSize', 15)
%         exportgraphics(fig, [(Cl{:}), '_', (Tx{:}), '_ORRHist.png'])
    end
end
%%
cntrs = (orrEdges(2:end) + orrEdges(1:end-1))/2;
for Cl= {'UMSCC22B','UMSCC47'}
    for Tx = {'NT', 'XT'}
            fig = figure; hold on
            fig.Position = [1 41 1600 783];
            title([Cl{:}, ' ', Tx{:}], 'FontSize',30)
            xlabel('ORR', 'FontSize',20)
            ylabel('Pixel Fraction', 'FontSize',20)
        for Ht = {'Baseline',  '24hpt', '48hpt'}
            try
                N = histcounts(histData.(Cl{:}).(Tx{:}).(matlab.lang.makeValidName(Ht{:})).ORRRaw, 25);
                plot(cntrs, N/sum(N), 'LineWidth', 4)           
            catch
                N =  histcounts(histData.(Cl{:}).NT.(matlab.lang.makeValidName(Ht{:})).ORRRaw, 25);
                plot(cntrs, N/sum(N), 'LineWidth', 4)
                continue
            end
        end
%         ylim([0, 0.05])
        legend({'Baseline',  '24hpt', '48hpt'}, 'FontSize', 15)
%         exportgraphics(fig, [(Cl{:}), '_', (Tx{:}), '_ORRHist.png'])
    end
end

%%
G = 0:1/99:1; S = 0:0.5/99:0.5;
for Cl = {'UMSCC22B','UMSCC47'}
    for Tx = {'XT'}
        for Ht = {'Baseline', '24hpt', '48hpt'}
            fig = figure; hold on
            fig.Position = [1 41 1600 783];
            title([Cl{:}, ' ', Tx{:}, ' ', Ht{:}], 'FontSize',30)
            xlabel('G', 'FontSize',20)
            ylabel('S', 'FontSize',20)
            try
                plt = surf(G, S, histData.(Cl{:}).(Tx{:}).(matlab.lang.makeValidName(Ht{:})).Phas);         
            catch
                plt = surf(G, S, histData.(Cl{:}).NT.(matlab.lang.makeValidName(Ht{:})).Phas);
                title([Cl{:}, ' ', Ht{:}], 'FontSize',30)
            end
            view(2)
            plt.LineStyle = 'None';
            shading interp
            cmap=colormap(jet);
            cmap(1,:) = 1;
            colormap(cmap);
            ax = gca;
            ax.YDir = 'normal'; ax.XLabel.String = 'G'; ax.YLabel.String= 'S';
            cnt = [0.5, 0];
            r = 0.5;
            x = 0:0.001:1;
            y = sqrt(r^2 - (x-cnt(1)).^2)+cnt(2);
            plot(x,y,'k', 'LineWidth', 5)
            cbar = colorbar; 
            cbar.Label.String = 'Pixel Fraction';
    
            exportgraphics(fig, [(Cl{:}), '_', (Tx{:}), '_', Ht{:}, '_Phasor.png'])
        end

    end
end
            
