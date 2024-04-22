function [G, S, photons] = getPhasorCoords(decay, params, varargin)

%% Parse inputs
if ~exist('decay', 'var') || isempty(decay)
    [decay_file, decay_path] = uigetfile( {'*.sdt'; '*.tif'}, 'Select decay file...');
    decay = [decay_path filesep decay_file];
end

% Input validators
isPars = @(x) isstruct(x) && all(isfield(x, {'dt', 'T', 'w'}));
validStyle = @(x) strcmpi(x, 'map') || strcmpi(x, 'mean');
validNumber = @(x) isnumeric(x) && isscalar(x) && (x>=0) && rem(x, 1)==0;
validFilter = @(x) all([size(x) == [1 2], arrayfun(validNumber, x)]) || validNumber(x);
isswitch= @(x) any(strcmpi(x, {'on', 'off'}));
validDecay = @(x) ((ischar(x) || isstring(x)) && isfile(x)) || (isnumeric(x) && ~ismatrix(x) && ~isscalar(x));
p = inputParser;
addRequired(p, 'decay', validDecay)
addRequired(p, 'params', isPars)
addOptional(p, 'calibrationwarning', 'on', isswitch)
addParameter(p, 'calibrationstyle', 'mean', validStyle);
addParameter(p, 'triggerdelay', false, @islogical)
addParameter(p, 'framedelay', 0, validNumber);
addParameter(p, 'timebinwidth', 1, validNumber);
addParameter(p, 'spatialbin', 1, validFilter);
addParameter(p, 'decaychannel', 1, validNumber);
addParameter(p, 'peakatzero', false, @isogical)
parse(p, decay, params, varargin{:})

%% Load Decay 
decay = p.Results.decay;
if ischar(decay) || isstring(decay) 
    [~, ~, ext] = fileparts(decay);
    switch ext
        case '.sdt'
            try
                decay = squeeze(double(SDTtoTimeStack(p.Results.decay, ...
                                        'framedelay', p.Results.framedelay, ...
                                      'triggerdelay', p.Results.triggerdelay, ...
                                           'channel', p.Results.decaychannel, ...
                                        'peakatzero', p.Results.peakatzero)));
            catch ME
                error('Exception occured while loading image. \n\n %s', getReport(ME))
            end
        case {'.tiff', '.tif'}
            try
                decay = double(tiffreadVolume(decay, 'PixelRegion', {[1 inf] [1 inf], [p.Results.decaychannel p.Results.decaychannel]}));
            catch
                done = false;
                ii = 1;
                while ~done
                    try
                        tempdecay(:,:,ii) = double(imread(decay, ii));
                        ii = ii+1;
                    catch ME
                        done = true;
                    end
                end
                decay = tempdecay(:,:,p.Results.decaychannel);
            end
    end
else
    if ndims(decay)>3
        decay = reshape(double(decay(:,:,p.Results.decaychannel,:)), size(decay, [1 2 ndims(decay)]));
    elseif ndims(decay) == 3 && p.Results.decaychannel ~=1
        error('Amiguous dimensions for channel selection \n. If third dimension is time, select channel 1. \n')
    end
end

%% Preprocess decay
% Spatial binning (defaults to no binning)
decay = imfilter(decay, ones(p.Results.spatialbin));

% Temporal binning (defaults to no binning)
decay = convn(decay, ones([1 1 p.Results.timebinwidth]));
decay = decay(:,:,p.Results.timebinwidth:p.Results.timebinwidth:end-mod(p.Results.params.T, p.Results.timebinwidth)); % This effectively performs a stride euqal to the bin width and ignores any final layers that did not have a full bin


%% Calculate Phasor coordinates
t = ones(size(decay)).*shiftdim(p.Results.params.dt*(1/2:1:p.Results.params.T-1/2), -1); % Create array of values of t for easy matops with peak at t=0.

% Calcuate phasor coords
g = sum(decay.*cos(p.Results.params.w*t), 3);
s = sum(decay.*sin(p.Results.params.w*t), 3);
photons = sum(decay, 3);

% Normalize Coords
g = g./photons;
s = s./photons;

% Adjust to calibration
style = lower(p.Results.calibrationstyle);
switch lower(p.Results.calibrationwarning)
    case 'off'
        G=g;
        S=s;
        return
    otherwise
        if ~isfield(p.Results.params, 'calibration') || ~isfield(p.Results.params.calibration, style)
             warning('Empty calibration structure. Raw G and S coordinates returned')
            return
        end
end
phi_adj = p.Results.params.calibration.(style)(p.Results.params.ch == p.Results.decaychannel).Phi;
m_adj = p.Results.params.calibration.(style)(p.Results.params.ch == p.Results.decaychannel).M;

[phi, m] = cart2pol(g, s);
[G, S] = pol2cart(phi+phi_adj, m.*m_adj);
end