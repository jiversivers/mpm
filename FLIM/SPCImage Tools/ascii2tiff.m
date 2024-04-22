function ascii2tiff(path)

% Collect all ascii files
files = dir([path, filesep, '*.asc']);

% Go through each file
fprintf(['Unpacking ',int2str(size(files,1)),' files..']);
for i = 1:size(files,1)
    try
        tmp = load([files(i).folder, filesep ,files(i).name],'-ascii');
        name = strsplit(files(i).name,'.asc');
        doublewrite(tmp,[files(i).folder,filesep ,name{1},'.tiff']);
    catch
        warning(['File ' files(i).name ' not read. Unexpected format.'])
    end
end
fprintf('Done\n')

%--- Internal function, ignore me ---%
function doublewrite(varargin)
% DOUBLEWRITE Write image with single precision to graphics file
%   DOUBLEWRITE(A,filename) writes image data A as a single precision
%   IEEEFP-formatted image to the file specified by filename. To that end, A
%   must be a single precision image.
%
%   DOUBLEWRITE(___,Name,Value) specifies additional parameters for ONLY tiff
%   image files, using one or more input arguments in the previous syntax.
%
% Created by: Alan Woessner (aewoessn@gmail.com)
% Last edited: Alan Woessner (3/19/2021)
%
% Maintained by: Kyle Quinn (kpquinn@uark.edu)
%                Quantitative Tissue Diagnostics Laboratory (Quinn Lab)
%                University of Arkansas   
%                Fayetteville, AR 72701

% Argument parser
input = argumentParser(varargin{:});

% Initialize tiff writer
if strcmp(input.writemode,'append')
    writer = Tiff(input.filename,'a');
else       
    writer = Tiff(input.filename,'w');
end
writer = setTags(writer,input.image);

% Write image
write(writer,input.image);
close(writer);

%--- Supporting functions ---%

%--- Argument Parser ---%
    function out = argumentParser(varargin)
        % Check to see if at least 2 inputs exist
        if (nargin < 2)
            error(message('MATLAB:imagesci:validate:wrongNumberOfInputs'));
        end
        
        % Initialize structure
        out = struct('image',[],'filename',[],'writemode','''');

        out.image = single(varargin{1});
        out.filename = varargin{2};
        
        % Parse an arbitrary amount of inputs
        for inpCount = 3:2:nargin
            out.(lower(varargin{inpCount})) = varargin{inpCount+1};
        end
        
        % Check to see if .tiff file extension is used
        split = strsplit(out.filename,'.');
        
        if ~contains(lower(split{end}),'tif')
            error(['Incorrect file extension used: .',split{end},', should be tiff extension']);
        end
    end

%--- TIFF Tag Setter ---%
    function writer = setTags(writer,image)
        % Set tags based on image
                
        % These depend on an RGB or grayscale image (currently, only
        % grayscale and RGB images are supported)
        if size(image,3) == 1
            % Grayscale image
            setTag(writer,'Photometric',Tiff.Photometric.MinIsBlack);
            setTag(writer,'BitsPerSample',32); % 32-bit
            setTag(writer,'SamplesPerPixel',1);
        elseif size(image,3) == 3
            % RGB image
            setTag(writer,'Photometric',Tiff.Photometric.RGB);
            setTag(writer,'BitsPerSample',32); % 32-bit
            setTag(writer,'SamplesPerPixel',3);
        else
            % Anything else
            error('Incorrect image size (should be either grayscale or RGB');
        end
        
        % These things are fixed
        setTag(writer,'ImageWidth',size(image,2));
        setTag(writer,'ImageLength',size(image,1));
        setTag(writer,'XResolution',1);
        setTag(writer,'YResolution',1);

        setTag(writer,'Compression',Tiff.Compression.None);
        setTag(writer,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
        setTag(writer,'ResolutionUnit',Tiff.ResolutionUnit.None);
        setTag(writer,'SampleFormat',Tiff.SampleFormat.IEEEFP);
    end
end
end