function [mask, note] = chooseMask(img, masks)
f1 = figure('WindowState', 'maximized');
img = double(img);
masks = double(masks);
mask = any(masks, 3);
note = 0;
for ii = 1:size(masks, 3)
    b = 0;
    while ~exist('b', 'var') || ~isempty(b)
        imagesc(mask.*img)
        [~, ~, b] = ginput(1);

        % Enter press -- Accept selection
        if isempty(b)
            continue

        % Space press -- Invert the mask
        elseif b == 32
            mask(logical(masks(:,:,ii))) = ~mask(logical(masks(:,:,ii)));

        % F/N Press -- Noisy, needs filtering (continue choosing)
        elseif b == 102 || b ==110
            note = 'F';

        % X/Esc Press -- Bad image/mask (will return the fcn)
        elseif b==27 || b==120
            mask = ones(size(img, [1 2]));
            note = 'X';
            warning('Mask or image marked as bad!')
            close(f1)
            return
        end
    end
end

close(f1)