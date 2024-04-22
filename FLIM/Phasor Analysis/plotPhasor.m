function varargout = plotPhasor(varargin)
%{
VARAGIN
    1: Axes to plot into; default is current ases
    2: G coordinates/array
    3: S coordinates/array
    Opt1: Edges of bins
    Opt2: Photon counts array (in same shape as G & S)
    ... Parameters for plot

VARAGOUT
    1: Phasor surface
    2: Histogram fractional counts array
    3: Bin centers
%}

p = inputParser;
edgy = {[-0.5, 1.5], [-0.25 0.75]};
isedgy = @(x) (size(x)==[1 2] || size(x)==[2 1]) && all(cellfun(@(y) all(size(y)==[1 2]) || all(size(y)==[2 1]), x));
isnotname = @(x) ismatrix(x) && isnumeric(x);
isaxis = @(x) contains(class(x), 'axis.Axes');
if ~isaxis(varargin{1})
    varargin = [{gca}, varargin(:)'];
end
addRequired(p, 'ax', isaxis)
addRequired(p, 'G', isnotname);
addRequired(p, 'S', isnotname);
addOptional(p, 'edges', edgy, isedgy)
addOptional(p, 'photons', ones(size(varargin{2})), isnotname)
addParameter(p, 'nbins', [128 128], isnotname);
addParameter(p, 'truncate', 1, isnotname);
addParameter(p, 'removeedge', false, @islogical);
addParameter(p, 'medfiltsize', [3 3], isnotname)
addParameter(p, 'medfiltcount', 0, @(x) isscalar(x) & isnotname(x))
addParameter(p, 'intensitythreshold', 0, isnotname)
parse(p, varargin{:});
G = p.Results.G;
S = p.Results.S;

% Drop origin vals ("infinite" lifetimes)
G(G==0 & S==0) = NaN;
S(G==0 & S==0) = NaN;

% Filter maps
if numel(p.Results.medfiltsize) == 1
    mf = [p.Results.medfiltsize p.Results.medfiltsize];
else
    mf = p.Results.medfiltsize;
end

if all([size(G), size(S)] ~= 1)
    for ii = 1:p.Results.medfiltcount
        G = medfilt2(G, mf);
        S = medfilt2(S, mf);
    end
elseif p.Results.medfiltcount > 0
    warning('Median Filter is not applied to vector inputs. Reshape G and S before input to use median filtering.')
end

% Flatten and threshhold (by default, 0 photon counts are removed)
G = G(p.Results.photons>p.Results.intensitythreshold);
S = S(p.Results.photons>p.Results.intensitythreshold);
G = reshape(G, numel(G), []);
S = reshape(S, size(G));

% Remove edge vals
if p.Results.removeedge
    edge = ~(G>-0.125 & G<1.125 & S>0 & S<0.75);
    G = G(~edge);
    S = S(~edge);
end

% Prep histogram bins
if numel(p.Results.nbins) == 1
    nbins = [p.Results.nbins p.Results.nbins];
else
    nbins = p.Results.nbins;
end
spaceG = 1/nbins(1);
ctrsG = p.Results.edges{1}(1)+0.5*spaceG:spaceG:p.Results.edges{1}(2)-0.5*spaceG;
spaceS = 0.5/nbins(2);   
ctrsS = p.Results.edges{2}(1)+0.5*spaceS:spaceS:p.Results.edges{2}(2)-0.5*spaceS;
varargout{3} = {ctrsG, ctrsS};

% Create 2D histogram surface
N = hist3([G, S], 'Ctrs', varargout{3});
N = 100*N/sum(N, 'all');
varargout{2} = N';

% Truncate high bins
N(N>p.Results.truncate) = p.Results.truncate;

% Plot
varargout{1} = surf(p.Results.ax, varargout{3}{:}, varargout{2});