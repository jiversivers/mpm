function phasorView(srf, cmap, refpts, params, circopt, flatten, logcon)

if ~exist('srf', 'var') || isempty(srf)
    srf = gca;
end

if exist('params', 'var')
    if isstruct(params)
        w = params.w;
    else 
        w = params;
    end
end

% Parse colormap input (or set to default)
if ~exist('cmap', 'var') || isempty(cmap)
    cmap= (colormap('jet')).*colormap('gray');
elseif ischar(cmap) || isstring(cmap)
    cmap = colormap(cmap);
end

view(2)
pbaspect([4 3 1])
if strcmp(srf.Type, 'axes')
    ax = srf;
    phasor = srf.Children(arrayfun(@(x) contains(class(x), 'Surface'), srf.Children));
else
    ax = srf.Parent;
    phasor = srf;
end
arrayfun(@(x) set(x, 'LineStyle', 'None'), phasor);
shading interp
colormap(cmap);

if exist('flatten', 'var') && (any(strcmp(flatten, {'flatten', 'on', 'true', 'flat'})) || flatten)
    phasor.CDataMode = 'manual';
    phasor.ZData = zeros(size(phasor.ZData));
end

if exist('logcon', 'var') && (any(strcmp(logcon, {'on', 'log', 'amp', 'contrast', 'logcon'})) || logcon)
    phasor.CData = log10(phasor.CData);
end

if ~exist("circopt", "var") || ~any(strcmp(circopt, {'circoff', 'off', 'nocircle', 'circleoff'}))
    z = arrayfun(@(x) get(x, 'Zdata'), phasor, 'UniformOutput', false);
    z = cat(3, z{:});
    circ = viscircles(ax, [0.5 0], 0.5);
    circ.Children(1).ZData = max(z, [], 'all')*ones(size(circ.Children(1).XData));
    circ.Children(2).ZData = circ.Children(1).ZData;
    if all(round(cmap(1,:))==0)
        circ.Children(1).Color = [0 0.25 0.75];
    else
        circ.Children(1).Color = [1 1 1];
    end
    circ.Children(2).Color = ~logical(round(cmap(1,:)));
end

if exist('refpts', 'var') && ~isempty(refpts) && exist('w', 'var')
    h = ishold(ax);
    hold(ax, 'on');
    [g, s] = lifetimetophasor(refpts*10^-9, params.w);
    scatter3(g, s, repmat(max(z+1, [], 'all'), 1, numel(refpts)), 15, [0.5 0 0], 'filled', 'o')

    [phi, m] = cart2pol(g-0.5, s);
    [g, s] = pol2cart(phi, m+0.03);
    lab = arrayfun(@(x) sprintf('%.2g ns', x), refpts, 'UniformOutput', false);
    text(g+0.5, s, repmat(max(z, [], 'all'), 1, numel(refpts)), lab, 'FontSize', 12, 'Color', [0.75 0 0], 'HorizontalAlignment', 'center')
    if ~h
        hold(ax, 'off')
    end
end
ax.XLabel.String = 'G';
ax.YLabel.String = 'S';
ax.YLim = [0 0.75];
ax.XLim = [0 1];
ax.XTick = [0 1];
ax.YTick = [0 0.5];