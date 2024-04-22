function [coloredMap, cmap] = createColorFitMap(fitParameter, photonCounts, edges)

% Create colormap
cmapLength = length(edges);
cmapJet = jet(100);
cmapJet = cmapJet(13:88,:);
cmap = zeros(cmapLength, 3);
for i = 1:3
    cmap(:,i) = interp1((1:length(cmapJet)),cmapJet(:,i),linspace(1,length(cmapJet),cmapLength));
end

% Map color onto parameter
[~,~,idx] = histcounts(fitParameter, edges);
idx(fitParameter>edges(end)) = cmapLength-1;
idx(fitParameter<edges(1)) = 0;
idx = idx+1;
coloredMap = ind2rgb(idx, cmap);

% Scale brightness by photon counts
allVals = nonzeros(sort(photonCounts(:)));
bot = allVals(round(0.05*numel(allVals)));
top = allVals(round(0.95*numel(allVals)));
photonCounts = (photonCounts-bot)./(top-bot);
photonCounts(photonCounts<0) = 0;
photonCounts(photonCounts>1) = 1;
coloredMap = coloredMap.*photonCounts;

% Convert to UINT8
coloredMap = uint8(round(coloredMap*255));