sim = readmatrix('ddAICsimResutls.xlsx');
cW = {'r', 'g', 'b', 'c', 'm'};
fig = figure;
fig.Position = [1921 221 1600 783];
for ii = 1:5
    pltData = [sim(:,4), sim(:,1), sim(:,11)];
    pltData = pltData(sim(:,10) == ii,:);

    scatter3(pltData(:,1), pltData(:,2), pltData(:,3), cW{ii})
    title(['Fits for ', num2str(ii), ' simulated components using ddAIC'])
    ylabel('Simulated Average of Averages')
    xlabel('Simulated Average of Standard Deviations')
    zlabel('Predicted Count of Subpopulations')
    ax = gca;
    ax.Position = [0.1300 0.1434 0.7322 0.7816];
    ax.CameraPosition = [-7, -3, 13];
    ax.XLim = [0 1];
    ax.YLim = [0 1];
    ax.ZLim = [1 5];

    % Save figure
    exportgraphics(fig, [num2str(ii), 'comp_ddAIC_figure.png'])

end