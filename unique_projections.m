folderPath = fullfile('Save');
fileNamePattern = '*_viablePoints_*.mat';

fileNames = dir(fullfile(folderPath, fileNamePattern));
fileNames = {fileNames.name};
nMod = numel(fileNames);
for iPar=nMod:-1:1
    fn = fileNames{iPar};
    loaded = load(fullfile(folderPath,fn));
    viablePointsArray(iPar) = loaded.viablePoints;
    % for duplicates removal
    projectionArray(iPar,:) = loaded.viablePoints.projection';
end
rmDuplicateProjections=true;  % Out of duplicates, take viable points with the largest sample
% Remove duplicates
if rmDuplicateProjections
    % Save init set for reference
    fileNames0 = fileNames;
    viablePointsArray0 = viablePointsArray;

    [uniqueProjectionArray, uniqueProjectionRowIdxs] = unique(projectionArray,'rows','stable');
    repeatsProjectionRowIdxs = setdiff(1:nMod, uniqueProjectionRowIdxs);
    repeatsProjectionArray = projectionArray(repeatsProjectionRowIdxs,:);
    nMod = numel(uniqueProjectionRowIdxs);
    for iPar=nMod:-1:1
        projectionIdx = uniqueProjectionRowIdxs(iPar);
        projection = projectionArray(projectionIdx,:);
        repeatsIdxs = find(ismember(repeatsProjectionArray, projection, 'rows'));
        %disp(repeatsProjectionRowIdxs(repeatsIdxs));
        maxSampSize = numel(viablePointsArray(projectionIdx).cost);
        for iRep=1:numel(repeatsIdxs)
            repeatIdx = repeatsProjectionRowIdxs(repeatsIdxs(iRep));
            repeatSampSize = numel(viablePointsArray(repeatIdx).cost);
            if repeatSampSize > maxSampSize
                maxSampSize = repeatSampSize;
                projectionIdx = repeatIdx;
            end
        end
        uniqueProjectionRowIdxs(iPar) = projectionIdx;
    end

    viablePointsArray = viablePointsArray(uniqueProjectionRowIdxs);
    fileNames = fileNames(uniqueProjectionRowIdxs);
end

    