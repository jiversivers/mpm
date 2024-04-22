function line = addMultiLifetimeLine(varargin)

% Parse positonal inputs
if ~contains(class(varargin{1}), 'matlab.graphics')
    phasor = gca;
    tau = varargin{1};
    w = varargin{2};
    varargin = varargin(3:end);
elseif contains(class(varargin{1}), 'Surface')
    phasor = varagin{1}.Parent;
    tau = varagin{2};
    w = varargin{3};
    varargin = varargin(4:end);
else
    phasor = varargin{1};
    tau = varargin{2};
    w = varargin{3};
    varargin = varargin(4:end);
end

h = ishold(phasor);
hold(phasor, 'on');

% Get ZData for line to be above phasor (visible over peak)
z = phasor.Children(contains(arrayfun(@class, phasor.Children, 'UniformOutput', false), 'Surface')).ZData;
z = max(z, [], 'all');

[g, s] = lifetimetophasor(tau, w);
line = plot(phasor, g', s', varargin{:});
for ii = 1:numel(line)
    line(ii).ZData = z*ones(size(line(ii).XData));
end

if ~h
    hold(phasor, 'off')
end