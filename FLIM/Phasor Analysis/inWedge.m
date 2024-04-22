function [tf, a, g, s] = inWedge(varargin)

% This funciton will return the TRUE for coordinates that are within a wedge
% of the universal circle from a set of test coordinates input. The first
% input defines the wedge as either a 1x3 vector of lifeitmes on the
% universal circle, with the first lifetime being the "pointy end" of the
% wedge and the otehr two being the range to fan out to, or as a 3x2 array
% of g and s coordinates, where the first row defines the pointy end and
% the last two rows define the points to fan out to. As an additinoal
% ouput, an array with the percent distance along the wedge (average of the
% distance along each edge of the wedge divided by the total length of the
% edge from the point to the circle) can be returned for each coordinate
% that falls into the wedge. This is analgous to the long-lifetime percent.
% The list of coordinates that are within the wedge from the test set can
% also be returned, if desired.

% Parse inputs
p = inputParser;
isWedge = @(x) all(size(x) == [1 3]) || all(size(x) == [3 2]);
addRequired(p, 'wedge', isWedge)
addRequired(p, 'coords', @isnumeric)
addOptional(p, 's', [], @isnumeric)
addParameter(p, 'w', [], @isnumeric)
parse(p, varargin{:});

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

if numel(p.Results.wedge) == 3
    if isempty(p.Results.w) 
        error('Angular frequency must be incldued when lifetimes are given as wedge points.')
    end
    [Gref, Sref] = lifetimetophasor(p.Results.wedge, p.Results.w);
else
    Gref = p.Results.wedge(:, 1);
    Sref = p.Results.wedge(:, 2);
end

% Shift coordinates so origin is at pointy end
x = G-Gref(1);
y = S-Sref(1);
gref = Gref-Gref(1);
sref = Sref-Sref(1);

% Convert new coords to polar to define wedge
[th, r] = cart2pol(x, y);
[ph, m] = cart2pol(gref, sref);

% Final check of values
tf = isWithin(th, ph(2:3), '=');
g = G.*tf;
s = S.*tf;
a = mean(cat(ndims(r)+1, (r.*tf)/m(2), (r.*tf)/m(3)), ndims(r)+1);


