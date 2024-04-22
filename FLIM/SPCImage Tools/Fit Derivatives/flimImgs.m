function flimData = flimImgs(flimData, limits)
% Makes intensity-weighted flim figures with upper and lower limits.
% flimData is a struct in the form of loadFlimFitResults outputs. limits is
% a 2x2 array with limit values. Row 1 is a1:a2 lower and upper limits, and
% row 2 is tauM lower and upper limits. If no limits are provided, default
% limits will be used. 0-5 for a1:a2 and 350-3000 for tauM.

% Default limits
if ~exist("limits", 'var')
    limits = [0 5 350 3000];
end

% Double all flimData
fields = fieldnames(flimData);
for ii = 1:numel(fields)
    flimData.(fields{ii}) = double(flimData.(fields{ii}));
end

% Calculate outputs (A1:A2 and TauM) maps
flimData.ratio = flimData.a1./flimData.a2;
flimData.tm = flimData.a1_per/100.*flimData.t1 + flimData.a2_per/100.*flimData.t2;

% Cutoff both maps to respective limits
ratio = flimData.ratio.*(flimData.ratio > limits(1) & flimData.ratio < limits(2));
tauM = flimData.tm.*(flimData.tm > limits(3) & flimData.tm < limits(4));

% Normalize maps and intensity image
int = flimData.photons/max(flimData.photons, [], 'all');

% Stretch maps to 8 bit
ratio = round(255 * ratio/5);
tauM = round(255 * tauM/5000);
int8 = round(255 * int);

% Create color map
cm = jet(340);
cm = cm(64:64+255,:);

% Map colors onto maps with intensity scaling
ratImg = zeros([size(ratio),3]);
tauImg = zeros([size(tauM),3]);
intImg = zeros([size(int), 3]);
for ii = 0:255
    for rgb = 1:3
        ratImg(:,:,rgb) = ratImg(:,:,rgb) + (cm(ii+1, rgb) * ((ratio == ii) .* int));
        tauImg(:,:,rgb) = tauImg(:,:,rgb) + (cm(ii+1, rgb) * ((tauM  == ii) .* int));
        intImg(:,:,rgb) = intImg(:,:,rgb) + (cm(ii+1, rgb) * ((int8  == ii) .* int));
    end
end

flimData.tmImg = tauImg;
flimData.tmMap = tauM;
flimData.ratioImg = ratImg;
flimData.ratioMap = ratio;