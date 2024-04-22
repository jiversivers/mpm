function trace = drawMask(img)

trace = {};
if size(img, 3) > 2
    image(gca, img(:,:,1:3))
else
    image(gca, mean(img, 3));
end

% Trace
KEY_PRESSED = 0;
ii = 0;
while ~KEY_PRESSED
    ii = ii+1;
    [~, x, y] = freehanddraw;
    trace{1, ii} = [x'; y'];   % Each section of mask goes into the next cell
    KEY_PRESSED = waitforbuttonpress;
end