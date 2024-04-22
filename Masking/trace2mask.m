function mask = trace2mask(trace, maskSize)
% Trace is an array or cell array of arrays of 2-row vectors with [x;y]
% line data to be converted to a logical polygon of size maskSize. If trace
% is a 2D cell array, each row is associated with a new image that will
% correspond to maskSize at that index. If maskSize is only size = 1, then
% the same maskSize will be used for all (in the case when all masked
% images in the batch are of the same size)

if ~isa(trace, 'cell')
    trace = {trace};
end

[xt, yt] = meshgrid(1:maskSize(2), 1:maskSize(1));
mask = zeros([maskSize, numel(trace)]);
for ii = 1:numel(trace)
    % Convert traces to masks
    stat = inpoly2([xt(:), yt(:)], trace{ii}');
    mask(:,:,ii) = reshape(stat, maskSize);
end
mask = logical(mask);