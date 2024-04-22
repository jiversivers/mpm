function [gp, sp] = projectPhasorToFit(varargin)

[G, S] = varargin{1:2};

% f is a linear fit object, m and b are sloep and y-int respectively
if nargin == 3
    f = varargin{3};
    m = f.p1;
    b = f.p2;
elseif nargin == 4
    m = varargin{3};
    b = varargin{4};
end

% Project
gp = (S+(G/m)-b)/(m+(1/m));
sp = m*gp+b;