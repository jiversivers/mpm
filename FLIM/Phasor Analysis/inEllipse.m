function [tf, g, s, roi] = inEllipse(varargin)

% Parse inputs
p = inputParser;
addRequired(p, 'ellipse')
addRequired(p, 'coords', @isnumeric)
addOptional(p, 's', [], @isnumeric)
parse(p, varargin{:});

% Parse coords
if ~isempty(p.Results.s)
    G = p.Results.coords;
    S = p.Results.s;
else 
    dimIdx = cell(size(size(p.Results.coords)));
    dimIdx(size(p.Results.coords) == 2) = {1};
    dimIdx(size(p.Results.coords) ~= 2) = {':'};
    G = p.Results.coords(dimIdx{:});
    dimIdx(size(p.Results.coords) == 2) = {2};
    S = p.Results.coords(dimIdx{:});
end

if ~all(size(G)==size(S))
    error('Incompatible G and S array sizes')
end

% If no ellipse was given, or draw was sleected, plot and draw ellipse
shape = size(G);
if isempty(p.Results.ellipse) || strcmp(p.Results.ellipse, 'draw')
    fig = figure;
    ax = axes(fig);
    % Prepare coords for plotting
    if any(shape == 1)
        f = factor(numel(G));
        G = reshape(G, prod(f(1:round(0.5*(numel(f))))), []);
        S = reshape(S, prod(f(1:round(0.5*(numel(f))))), []);
    end
    phasor = plotPhasor(ax, G, S, 'nbins', [128 128], 'medfiltcount', 1);
    phasorView(phasor, [], [], [], [], 'flatten')
    hold on
    ellipse = drawellipse;
    wait(ellipse);
    roi = ellipse.Vertices;
    close(fig);
    % Reshape back into original shape
    G = reshape(G, shape);
    S = reshape(S, shape);
elseif isa(p.Results.ellipse, 'matlab.graphics.axis.Axes') || isa(p.Results.ellipse, 'matlab.ui.Figure') || isa(p.Results.ellipse, 'matlab.graphics.chart.primitive.Surface')
    % Prepare coords for plotting
    if any(shape == 1)
        f = factor(numel(G));
        G = reshape(G, prod(f(1:round(0.5*(numel(f))))), []);
        S = reshape(S, prod(f(1:round(0.5*(numel(f))))), []);
    end
    hold on
    ellipse = drawellipse;
    wait(ellipse);
    roi = ellipse.Vertices;
    % Reshape back inot original shape
    G = reshape(G, shape);
    S = reshape(S, shape);

else
    roi = p.Results.ellipse;
end

% Final in/out checks
tf = inpoly2([G(:), S(:)], roi);
tf = reshape(tf, shape);
g = G.*tf;
s = S.*tf;
