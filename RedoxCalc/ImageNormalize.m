function normedImage = ImageNormalize(img, metadata, option)
if ~exist("option","var")
    option = 'direct';
end

if ~isstruct(metadata.Power)
    warning('Cannot normalize with missing power.')
    normedImage = NaN(size(img));
    return
end

switch lower(option)
    case {'f', 'fluorescein'}
        opts = struct();
        
        % Prepare metadata
        opts.RefAtten = cell2mat(keys(metadata.Power.Laser));
        opts.RefPow = values(metadata.Power.Objective);
        opts.RefAtten = opts.RefAtten(~cellfun(@ismissing, opts.RefPow));
        opts.RefPow = cell2mat(opts.RefPow(~cellfun(@ismissing, opts.RefPow)));
%         opts.Attenuation = opts.RefAtten(abs(metadata.Pockels-opts.RefAtten) == min(abs(metadata.Pockels-opts.RefAtten)));
        opts.Attenuation = metadata.Pockels;
        opts.dateTaken = metadata.Date;
        
        if strcmp(metadata.Laser, 'InsightX3')
            opts.laserFlag = 'Upright1';
        else
            opts.laserFlag = 'Inverted1';
        end
        
        normedImage = zeros(size(img));
        for CH = 1:size(img, 3)
            opts.PMTChannel = CH;
            opts.PMT = metadata.Gain(CH);
            % Normalize
            normedImage(:,:,CH) = fluorescein_normalization_master(img(:,:,CH), opts);
        end
    case {'d', 'direct'}
        % Prepare metadata
        att = metadata.Pockels;
        refAtt = cell2mat(keys(metadata.Power.Objective));
        refPwr = values(metadata.Power.Objective);
        refAtt = refAtt(~cellfun(@ismissing, refPwr));
        refPwr = cell2mat(refPwr(~cellfun(@ismissing, refPwr)));
        pmt = metadata.Gain;

        % Find the closest ref attenuation and pull that power measurement
        pwr = refPwr(abs(att-refAtt) == min(abs(att-refAtt)));

        % Normalize
        normedImage = double(img)./(shiftdim(pmt, -1)*pwr.^2);
    otherwise
        error('Select a valid normalization option or leave blank for direct normalization')
end
