%% Count subPops

% Import excel spreadsheet into cell array of tables
fn = 'ORR_freqCurve.xlsx';
xlSheets = sheetnames(fn);
histFits = readtable(fn, 'Sheet', xlSheets{1});

animalCol = histFits.AnimalID;
sampleCol = histFits.Sample_;
fovCol = histFits.FOV_;


%% Counting pops for each level of hierarchy

animalID = unique(animalCol);

% Total # of animals from the study
subPopStats.animalCount = length(animalID);

% Initiate global pop & fov counters
subPopStats.popCount = 0;
subPopStats.fovCount = 0;   

for ii = 1:length(animalID)
    % Isolate rows for animal
    animalRows = strcmp(animalCol, animalID{ii});
    sampleID = unique(sampleCol(animalRows));

    % Total # of samples from the animal
    subPopStats.(matlab.lang.makeValidName(animalID{ii})).samCount = length(sampleID);

    % Initiate animal pop & fov counters
    subPopStats.(matlab.lang.makeValidName(animalID{ii})).popCount = 0;
    subPopStats.(matlab.lang.makeValidName(animalID{ii})).fovCount = 0;

    for jj = 1:length(sampleID)
        % Isolate rows for this animal-sample
        sampleRows = animalRows & strcmp(sampleCol, sampleID{jj});
        fovID = unique(fovCol(sampleRows));

        % Total # of fovs in the sample
        subPopStats.(matlab.lang.makeValidName(animalID{ii})).(matlab.lang.makeValidName(sampleID{jj})).fovCount = length(fovID);

        % Intiate sample population counter
        subPopStats.(matlab.lang.makeValidName(animalID{ii})).(matlab.lang.makeValidName(sampleID{jj})).popCount = 0;
        
        for kk = 1:length(fovID)
            % Pops in current fov (Isolate rows for animal-sample-fov)
            subPopStats.(matlab.lang.makeValidName(animalID{ii})).(matlab.lang.makeValidName(sampleID{jj})).(matlab.lang.makeValidName(fovID{kk})).popCount = sum(sampleRows & strcmp(fovCol, fovID{kk}));

            % Increment pop counter for current sample by current fov pops
            subPopStats.(matlab.lang.makeValidName(animalID{ii})).(matlab.lang.makeValidName(sampleID{jj})).popCount = subPopStats.(matlab.lang.makeValidName(animalID{ii})).(matlab.lang.makeValidName(sampleID{jj})).popCount + subPopStats.(matlab.lang.makeValidName(animalID{ii})).(matlab.lang.makeValidName(sampleID{jj})).(matlab.lang.makeValidName(fovID{kk})).popCount;
        end
        
        % Increment pop count for current animal by current sample pops
        subPopStats.(matlab.lang.makeValidName(animalID{ii})).popCount = subPopStats.(matlab.lang.makeValidName(animalID{ii})).popCount + subPopStats.(matlab.lang.makeValidName(animalID{ii})).(matlab.lang.makeValidName(sampleID{jj})).popCount;
        
        % Increment fov count for current animal by sample fovs
        subPopStats.(matlab.lang.makeValidName(animalID{ii})).fovCount = subPopStats.(matlab.lang.makeValidName(animalID{ii})).fovCount + subPopStats.(matlab.lang.makeValidName(animalID{ii})).(matlab.lang.makeValidName(sampleID{jj})).fovCount;
    end
    
    % Increment pop count for current study by current animal pops
    subPopStats.popCount = subPopStats.popCount + subPopStats.(matlab.lang.makeValidName(animalID{ii})).popCount;
    
    % Increment fov count for current study by current animal pops
    subPopStats.fovCount = subPopStats.fovCount + subPopStats.(matlab.lang.makeValidName(animalID{ii})).fovCount;

end

%% Calculate stats at each level

clear sampleID fovID

subPopStats.meanPopPerAnimal = subPopStats.popCount / subPopStats.animalCount;
subPopStats.meanPopPerFov = subPopStats.popCount / subPopStats.fovCount;

gloVarNum = 0;
for ii = 1:length(animalID)
    subPopStats.(animalID{ii}).meanPopPerFOV = subPopStats.(animalID{ii}).popCount / subPopStats.(animalID{ii}).fovCount;
    sampleID = fieldnames(subPopStats.(animalID{ii}));
    aniVarNum = 0;
    for jj = 1:length(sampleID)
        if isstruct(subPopStats.(animalID{ii}).(sampleID{jj}))
            subPopStats.(animalID{ii}).(sampleID{jj}).meanPopPerFOV = subPopStats.(animalID{ii}).(sampleID{jj}).popCount / subPopStats.(animalID{ii}).(sampleID{jj}).fovCount;
            fovID = fieldnames(subPopStats.(animalID{ii}).(sampleID{jj}));
            samVarNum = 0;
            for kk = 1:length(fovID)
                if isstruct(subPopStats.(animalID{ii}).(sampleID{jj}).(fovID{kk}))
                    numInc = (subPopStats.(animalID{ii}).(sampleID{jj}).(fovID{kk}).popCount - subPopStats.(animalID{ii}).(sampleID{jj}).meanPopPerFOV)^2;
                    samVarNum = samVarNum + numInc;
                    aniVarNum = aniVarNum + numInc;
                    gloVarNum = gloVarNum + numInc;
                end
            end
        subPopStats.(animalID{ii}).(sampleID{jj}).varInSam = samVarNum/(subPopStats.(animalID{ii}).(sampleID{jj}).fovCount-1);
        end
    end
    subPopStats.(animalID{ii}).varInAnimal = aniVarNum/(subPopStats.(animalID{ii}).fovCount-1);
end

subPopStats.varInStudy = gloVarNum / subPopStats.fovCount;

%% Save

save('ORR_subPopStats_fromCurveStats', 'subPopStats')
writetable(NestedStruct2table(subPopStats),'ORR_subPopStats_fromCurveStats.xlsx')