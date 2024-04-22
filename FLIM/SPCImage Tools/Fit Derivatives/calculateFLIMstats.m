function flimStats = calculateFLIMstats(flimData)
%{
This function takes the fit data structure (from the loadFitData function
above) and adds mean and standard deviation statistics ('Avg_<parameter
name>' and 'StD_<parameter name>' respectively for dot-indexing) to the
structure. It also adds derivative parameters, namely, mean lifetime and
the free to bound ratio ('tm' and 'a1a2ratio' respectively for
dot-indexing).
%}

fname = fieldnames(flimData);

% % Calculate 2nd order data (tm and a1a2) ratio if not included in input
% if ~contains('tm', fname) || ~contains('a1a2ratio', fname)
%     flimData = flimCombination.flimData;
% end
    

% Calculate all stats for input data
flimStats = struct();
for ii = 1:length(fname)
    if ~(contains(fname{ii}, 'Img') || contains('image', fname{ii}) || contains('IRF', fname{ii}) || contains('curve', fname{ii}) || contains('trace', fname{ii}) || contains('legend', fname{ii}))
        flimStats.(['Avg_', fname{ii}]) = mean(flimData.(fname{ii}), 'all', 'omitnan');
        flimStats.(['StD_', fname{ii}]) = std(flimData.(fname{ii}), 0, 'all', 'omitnan');
    end
end

% Calculate bulk stats for second order data (tm and a1a2 ratio)
flimStats.Bulk_Avg_tm = flimStats.Avg_a1_per/100 * flimStats.Avg_t1 + flimStats.Avg_a2_per/100 * flimStats.Avg_t2;
flimStats.Bulk_StD_tm = flimStats.Avg_a1_per/100 * flimStats.StD_t1 + flimStats.Avg_a2_per/100 * flimStats.StD_t2;
flimStats.Bulk_Avg_ratio = flimStats.Avg_a1/flimStats.Avg_a2;
