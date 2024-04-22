startDir = pwd;

% Find FLIM data folders
flimDirs = FolderFinder('LifetimeData_Cycle00001_000001.sdt');

parts = pathpieces(flimDirs);
%% 
outTable = cell([length(flimDirs), 2]);
for ii = 1:length(flimDirs)
    cd(flimDirs{ii})
    
    % Reading current lambda .XML file
    xml = xml2struct_nonnative(dir(['*.xml']).name);
    try
        if xml.PVScan.PVStateShard.PVStateValue{1, 10}.IndexedValue.Attributes.description == 'PockelsMaiTai'
            outTable(ii, :) = {[parts{ii,end-1}, '-', parts{ii,end}], 'Upright'};
        end 
    catch
        try
            if xml.PVScan.PVStateShard.PVStateValue{1, 11}.IndexedValue.Attributes.description == 'MaiTai'
                outTable(ii, :) = {[parts{ii,end-1}, '-', parts{ii,end}], 'Inverted'};
            end
        catch
            if  xml.PVScan.PVStateShard.PVStateValue{1, 8}.IndexedValue.Attributes.description == 'PockelsMaiTai'
                outTable(ii, :) = {[parts{ii,end-1}, '-', parts{ii,end}], 'Upright'};
            else
                outTable(ii, :) = {[parts{ii,end-1}, '-', parts{ii,end}], 'Unrecognized'};
            end
        end
    end
end

cd(startDir)
writecell(outTable, 'LaserUsed.xlsx')



