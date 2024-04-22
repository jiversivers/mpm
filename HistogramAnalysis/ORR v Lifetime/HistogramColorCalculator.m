function [C, n] = HistogramColorCalculator(Var1, Var2, colorBasis, binCtrs)

edgs = cellfun(@(x) [x-((x(2)-x(1))/2) x(end)+((x(2)-x(1))/2)], binCtrs, 'UniformOutput', false);
C = zeros(flip(cellfun(@numel, binCtrs)));
n = zeros(flip(cellfun(@numel, binCtrs)));
    for x = 1:numel(edgs{1})-1
        inX = Var1>=edgs{1}(x) & Var1<edgs{1}(x+1);
        for y = 1:numel(edgs{2})-1
            inY = Var2>=edgs{2}(y) & Var2<edgs{2}(y+1);
            c = colorBasis(inX & inY);
            C(y, x) = mean(c, 'omitnan');
            n(y, x) = numel(c);
        end
    end