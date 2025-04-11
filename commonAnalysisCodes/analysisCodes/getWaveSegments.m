function [allDir,allUniqueDirs,newBounds,speed,tempFreq,clusters,pgd,emptyCells] = getWaveSegments(outputsTW,timeVals,wobbleReq)
% set boundryLims as a vector of times outside of which direction vectors
% will not be considered
% if timeFlag = 1, append time values, else append the indices instead
boundryLims = [0.25 0.75];
boundryLims = [dsearchn(timeVals',boundryLims(1)),dsearchn(timeVals',boundryLims(2))];
%% find the travelling wave indices for all trials and plot
pgd = zeros(numel(outputsTW),length(timeVals));
directions = zeros(numel(outputsTW),length(timeVals));
tempFreq = zeros(numel(outputsTW),length(timeVals)-1);
clusters = zeros(numel(outputsTW),length(timeVals));
speed = zeros(numel(outputsTW),length(timeVals));

for i = 1:numel(outputsTW)
    pgd(i,:) = outputsTW{i}.pgd;
    directions(i,:) = outputsTW{i}.direction;
    speed(i,:) = outputsTW{i}.speed;
    tempFreq(i,:) = outputsTW{i}.tempFreq;
    clusters(i,:) = outputsTW{i}.clusters;
end
tempFreq = cat(2,zeros(size(tempFreq,1),1),tempFreq);

pgd(pgd<0) = 0;
pgd(isnan(pgd)) = 0;
directions(pgd==0) = nan;
directions(:,1:boundryLims(1)-1) = nan;
directions(:,boundryLims(2)+1:length(timeVals)) = nan;

allDir = nan(size(directions,1),length(timeVals));
newBounds = cell(1,size(directions,1));
allUniqueDirs = cell(1,size(directions,1));

%% diverge into type of segmentation required
if wobbleReq==1
    for i = 1:size(directions,1)
        boundries = [];
        directionInt = zeros(1,length(timeVals));
        directionInt(~isnan(directions(i,:))) = 1;
        directionData = find(directionInt==0);
        directionEpochs = find(diff(directionData)>1);
        boundries = cat(2,boundries,[directionData(directionEpochs)+1;directionData(directionEpochs+1)-1]); 
        boundries(:,diff(boundries)<5) = [];
        for k = 1:size(boundries,2)
            allDir(i,boundries(1,k):boundries(2,k)) = directions(i,boundries(1,k):boundries(2,k));
            allUniqueDirs{i,k} = circ_mean(unique(directions(i,boundries(1,k):boundries(2,k)))');
        end
        newBounds{i} = boundries;
    end
    emptyCells = [];
    speed(isnan(allDir)) = nan;
    tempFreq(isnan(allDir)) = nan;
    clusters(isnan(allDir)) = nan;
    pgd(isnan(allDir)) = nan;

elseif wobbleReq==2
    wobbleLim = 10;
    for i = 1:size(directions,1)
    boundries = [];
    directionInt = zeros(1,length(timeVals));
    directionInt(~isnan(directions(i,:))) = 1;
    directionData = find(directionInt==0);
    directionEpochs = find(diff(directionData)>1);
    boundries = cat(2,boundries,[directionData(directionEpochs)+1;directionData(directionEpochs+1)-1]); 
    boundries(:,diff(boundries)<40) = [];
    revisedBounds = [];
        for k = 1:size(boundries,2)
            tempBounds = directions(i,boundries(1,k):boundries(2,k));
            tempBoundsInd = boundries(1,k):boundries(2,k);
            tempDiff = abs(diff(tempBounds));
            if max(tempDiff)<deg2rad(wobbleLim)
                revisedBounds = cat(2,revisedBounds,[boundries(1,k);boundries(2,k)]);
            else
            stpPts = find(tempDiff>deg2rad(wobbleLim))+1;
            stpPts = [1 stpPts;stpPts-1 numel(tempBounds)];
            revisedBounds = cat(2,revisedBounds,[tempBoundsInd(stpPts(1,:));tempBoundsInd(stpPts(2,:))]);
            end
            revisedBounds(:,diff(revisedBounds)<40) = [];
        end 
        uniqueDirs = [];
        for j = 1:size(revisedBounds,2)
            allDir(i,revisedBounds(1,j):revisedBounds(2,j)) = directions(i,revisedBounds(1,j):revisedBounds(2,j));
            uniqueDirs = cat(2,uniqueDirs,circ_mean(directions(i,revisedBounds(1,j):revisedBounds(2,j))'));
        end
        newBounds{i} = revisedBounds;
        allUniqueDirs{i} = uniqueDirs;
    end
speed(isnan(allDir)) = nan;
tempFreq(isnan(allDir)) = nan;
clusters(isnan(allDir)) = nan;
pgd(isnan(allDir)) = nan;

emptyCells = find(cellfun(@isempty,newBounds));
newBounds(emptyCells) = [];
speed(emptyCells,:) = [];
tempFreq(emptyCells,:) = [];
clusters(emptyCells,:) = [];
pgd(emptyCells,:) = [];

else
    wobbleLim = 5;
    for i = 1:size(directions,1)
        boundries = getLims(directions(i,:),1);
        boundries(:,diff(boundries)<50) = [];
        revisedBounds = [];
        for k = 1:size(boundries,2)
            tempBounds = directions(i,boundries(1,k):boundries(2,k));
            boundaryIndices = boundries(1,k):boundries(2,k);
            waveIndices = getSegmentedWaves(tempBounds,boundaryIndices,wobbleLim);
            revisedBounds = cat(2,revisedBounds,waveIndices);
        end 
        revisedBounds(:,diff(revisedBounds)<50) = [];

        if ~isempty(revisedBounds)
        uniqueDirs = [];
            for j = 1:size(revisedBounds,2)
                allDir(i,revisedBounds(1,j):revisedBounds(2,j)) = directions(i,revisedBounds(1,j):revisedBounds(2,j));
                uniqueDirs = cat(2,uniqueDirs,circ_mean(directions(i,revisedBounds(1,j):revisedBounds(2,j))'));
            end
        newBounds{i} = revisedBounds;
        allUniqueDirs{i} = uniqueDirs;
        end
    end

    speed(isnan(allDir)) = nan;
    tempFreq(isnan(allDir)) = nan;
    clusters(isnan(allDir)) = nan;
    pgd(isnan(allDir)) = nan;
    emptyCells = find(cellfun(@isempty,newBounds));
    newBounds(emptyCells) = [];
    speed(emptyCells,:) = [];
    tempFreq(emptyCells,:) = [];
    clusters(emptyCells,:) = [];
    pgd(emptyCells,:) = [];
end
end

%% additional functions
%%
function [waveIndices] = getSegmentedWaves(data,boundaryIndices,devAngle)
dataFull = [nan data nan];
uniqueAngles = unique(dataFull);
uniqueAngles(isnan(uniqueAngles)) = [];

segWaveBounds = [];
inMat = zeros(length(segWaveBounds),length(segWaveBounds));
for i = 1:numel(uniqueAngles)
    data = dataFull;
    initAngle = uniqueAngles(i);
    deviation = deg2rad(devAngle);
    dataDiff = abs(data-initAngle);
    data(dataDiff>deviation) = 0;
    data(data==0) = nan;

    boundsInt = getLims(data,1);
    segWaveBounds = cat(2,segWaveBounds,boundsInt);
end

segWaveBounds = segWaveBounds-1;
[~,b] = sort(diff(segWaveBounds),'descend');
sortedBounds = segWaveBounds(:,b);
durations = diff(sortedBounds);

% finalBounds = zeros(length(sortedBounds),length(sortedBounds)-1);
if size(sortedBounds,2)>1
    for i = 1:length(sortedBounds)
        wave1 = sortedBounds(1,i):sortedBounds(2,i);
        for j = 1:length(sortedBounds)
            wave2 = sortedBounds(1,j):sortedBounds(2,j);
            inMat(i,j) = numel(intersect(wave1,wave2));
        end
    end
    inMat(find(diag(diag(inMat)))) = 1;
    inMat(inMat==0) = 0;

    for i = 1:length(sortedBounds)
        selection = inMat(i,:);
        overlapWaves = find(selection>0);
        if ~isempty(overlapWaves)
            for j = 1:numel(overlapWaves)
                if durations(i)== durations(overlapWaves(j))
                    inMat(i,overlapWaves(j)) = 1;
                    inMat(overlapWaves(j),i) = 1;
                elseif durations(i)>durations(overlapWaves(j))
                    inMat(overlapWaves(j),:) = nan;
                    inMat(:,overlapWaves(j)) = nan;
                end
            end
        end
    end
    nonOverlappingWaves = sortedBounds(:,find(diag(inMat)==1));
else 
    nonOverlappingWaves = sortedBounds;
end
waveIndices = zeros(2,size(nonOverlappingWaves,2));
for i = 1:size(nonOverlappingWaves,2)
    waveIndices(:,i) = [boundaryIndices(nonOverlappingWaves(1,i));boundaryIndices(nonOverlappingWaves(2,i))];
end
end

function [boundries] = getLims(data,flag)
if flag==1
    % for segmenting when required segments are sparse (intervening zeros)
    boundries = [];
    dataInt = zeros(1,length(data));
    dataInt(~isnan(data)) = 1;
    dataSig = find(dataInt==0);
    dataEpochs = find(diff(dataSig)>1);
    boundries = cat(2,boundries,[dataSig(dataEpochs)+1;dataSig(dataEpochs+1)-1]); 
else
    % for segmenting continuous segments 
    % dataDiff = diff(data);
    dataDiff = [1,abs(diff(data)),1];
    dataDiff(dataDiff~=0) = 1;
    boundries = find(dataDiff);
    boundries = [boundries(1:end-1);boundries(2:end)-1];
end
end