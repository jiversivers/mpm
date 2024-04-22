function t = labelUniversalCirlce(g,s,lab, varargin)
[th, r] = cart2pol(g-0.5,s);
[g, s] = pol2cart(th, r+0.05);
g = g+0.5;
t = text(g, s, lab, varargin{:}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
