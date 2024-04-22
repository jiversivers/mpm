% This script will automatically find and rename B&H files for batch
% processing. It will copy them into a seperate directory in the selected
% directory. It will only find files further down the directory hierarchy
% from the selected root. It will automatically copy and rename all the
% files based on the information in the directory names beneath the root.
% If the wavelength of acquisition was with excitation of 755 and/or 855 nm
% and  this is included in the image directory name, the files will be
% sorted into subfolders of these names (to prepared for individual
% wavelength-based batches). If this inofrmation is not evident in the
% image directory name or the images were acquired ata different
% exctiation, the files will be copied to a subdirectory titled "Other."

% Set up directory for data copy
startDir = uigetdir(pwd, 'Select root directory of data tree');
saveDir = uigetdir(startDir, 'Select save directory');
d = numel(strsplit(startDir, filesep))+1;

terminals = FolderFinder('LifetimeData_Cycle00001_000001.sdt', startDir);
np = pathpieces(terminals);
np([contains(np(:, end), '755')], end+1) = {'755'};
np([contains(np(:, end-1), '855')], end) = {'855'};
np([cellfun(@isempty, np(:,end))], end) = {'Other'};
subs = unique(np(:,end));
for ii = 1:numel(subs)
    if ~exist([saveDir filesep subs{ii}], 'dir')
        mkdir([saveDir filesep subs{ii}])
    end
end

for ii = 1:length(terminals)
    % Extract name parts from dir hierarchy below root and form into file names
    newName = namemaker(np(ii, d:end-1), '_');
    try
        if ~exist([saveDir np{ii, end} filesep newName,'.sdt'], 'file') 
            copyfile([terminals{ii} filesep 'LifetimeData_Cycle00001_000001.sdt'],[saveDir filesep np{ii, end} filesep newName,'.sdt'])
            disp(['Copying and renaming B&H for ' newName, '.'])
        else
            disp(['B&H File named ' newName,' already exists in FLIM_Batch. Copy skipped.'])
        end
    catch ME
        error(getReport(ME))
        disp(['B&H File named ' newName,' could not be sorted. Copy skipped.'])
    end
end

cd(saveDir)
cd ..\