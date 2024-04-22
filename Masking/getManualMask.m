function [totalMasks, masks] = getManualMask(impath)
%% 
% This fucntion will dislay the image and allow you to draw a mask in
% multiple sections (for donut-type masking) Double click to complete a
% section and single click to begin the next one. Depending on the image
% and mask size, the time between the drawing and mask selection step
% (described below) can be large, so to improve performance for batch
% processing, each step is run seperately for the entire batch. This will
% allow all masks to be drawn before selection. You will then choose
% whether it is a positive or negative mask using the space bar. Press
% enter to move onto the next section or complete selection (for the last
% selection).

% Batch or single?
if ~isa(impath, "cell")
    impath = {impath};      % To be treated like batch of n=1
end

% Draw mask(s)
masks = drawMask(vim);

% Invert/Revert masks
masks = cellfun(@(x,y) chooseMask(x, y), vim, masks, 'UniformOutput', false);

% Combine masks
totalMasks = cellfun(@(x) sum(x, 3), masks, 'UniformOutput', false);