function figs = makeFLIMFigs(flimData, limits)

% Make image arrays
flimData = flimImgs(flimData);

% Create A1:A2 figure
ratFig = figure;
ratFig.Position = [1921 221 783 783];
image(flimData.ratioImg)
axis image
colormap jet
c1 = colorbar;
caxis(gca, limits(1:2));
c1.Label.Interpreter = 'latex';
c1.Label.String = '$\frac{\alpha_1}{\alpha_2}$';
c1.Label.FontSize = 18;
c1.Label.Rotation = 0;
c1.Label.Position = [3 mean(limits(1:2)) 0];
axis off

% Create TauM figure
tauFig = figure;
tauFig.Position = [1921 221 783 783];
image(flimData.tmImg)
axis image
colormap jet
c2 = colorbar;
caxis(gca, limits(3:4));
c2.Label.Interpreter = 'latex';
c2.Label.String = '$\tau$ (ps)';
c2.Label.FontSize = 18;
c2.Label.Rotation = 0;
c2.Label.Position = [3.5 mean(limits(3:4)) 0];
axis off

% % Create Phasor figure
% phaFig = plotPhasor(flimData.G, flimData.S);

figs = struct('A1A2', ratFig, 'TauM', tauFig);