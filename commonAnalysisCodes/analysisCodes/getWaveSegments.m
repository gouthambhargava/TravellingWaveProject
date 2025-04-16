function [waveVector,uniqueDirs,waveBounds] = getWaveSegments(outputsTW,timeVals,wobbleLim,segOption,boundryLims, lengthLimit)
% set boundryLims as a vector of times outside of which direction vectors will not be considered
% Inputs: -outputsTW - output from getTWCircParams which contains the PGD, direction, speed, tempFreq and clusters for a single trial
%         -timeVals - vector of time values for the given trial
%         -wobbleLim - sets a threshold of an angle(in deg). Consecutive
%         values in direction under this threshold are considered part of a
%         single wave. Leave this option empty if segmentation according to
%         direction values is not required
%         -segOption - 
%         -boundryLims - give time limits, ex [0.25 0.75]. Waves will be identified only within this limit  
%         -lengthLimit - set a length value (in ms). Only waves above this
%         limit be considered. 

% define some parameters
boundryLims = dsearchn(timeVals',boundryLims(1)):dsearchn(timeVals',boundryLims(2)); %indices of the specified time limit
Fs = 1/(timeVals(2)-timeVals(1));
lengthLimit = lengthLimit*Fs*10^-3;
direction = outputsTW.direction;
direction(setdiff(1:length(timeVals),boundryLims)) = nan;
%% diverge into type of segmentation required

if segOption==1 % get wave segments where segments are purely seperated by significant PGD. Segmentation is very lenient 
    [waveVector,waveBounds] = simpleWaveSegments(direction,lengthLimit);
    uniqueDirs = [];
    for j = 1:size(waveBounds,2)
        uniqueDirs = cat(2,uniqueDirs,circ_mean(direction(waveBounds(1,j):waveBounds(2,j))'));
    end

elseif segOption==2 % Imposes an additional criteria on wave segmentation, based on wobbleLim.
    % waves within this earlier segmentation are further degraded based on
    % the angular difference between successive time points. Each time
    % point can vary only within the wobbleLim. Variation beyond this limit
    % marks the start of a new wave.
    

    [newDirection,boundries] = simpleWaveSegments(direction,lengthLimit);

    for k = 1:size(boundries,2)
        tempBounds = newDirection(boundries(1,k):boundries(2,k));
        tempBoundsInd = boundries(1,k):boundries(2,k);
        tempDiff = abs(diff(tempBounds));
        revisedBounds = [];
        if max(tempDiff)<deg2rad(wobbleLim)
            revisedBounds = cat(2,revisedBounds,[boundries(1,k);boundries(2,k)]);
        else
            stpPts = find(tempDiff>deg2rad(wobbleLim))+1;
            stpPts = [1 stpPts;stpPts-1 numel(tempBounds)];
            revisedBounds = cat(2,revisedBounds,[tempBoundsInd(stpPts(1,:));tempBoundsInd(stpPts(2,:))]);
        end
        revisedBounds(:,diff(revisedBounds)<lengthLimit) = [];
    end
    waveVector = nan(1,length(direction));
    uniqueDirs = [];
    for j = 1:size(revisedBounds,2)
        waveVector(revisedBounds(1,j):revisedBounds(2,j)) = direction(revisedBounds(1,j):revisedBounds(2,j));
        uniqueDirs = cat(2,uniqueDirs,circ_mean(direction(revisedBounds(1,j):revisedBounds(2,j))'));
    end
    waveBounds = revisedBounds;
    
else % ensures that wave segments from option 1 are further fragmented so that the total amount of variation between
    % the start and end of a wave is within the wobbleLim
    [~,boundries] = simpleWaveSegments(direction,lengthLimit);
    revisedBounds = [];
    for k = 1:size(boundries,2)
        tempBounds = direction(boundries(1,k):boundries(2,k));
        boundaryIndices = boundries(1,k):boundries(2,k);
        waveIndices = getSegmentedWaves(tempBounds,boundaryIndices,wobbleLim);
        revisedBounds = cat(2,revisedBounds,waveIndices);
    end
    revisedBounds(:,diff(revisedBounds)<lengthLimit) = [];
    
    waveVector = nan(1,length(direction));
    if ~isempty(revisedBounds)
        uniqueDirs = [];
        for j = 1:size(revisedBounds,2)
            waveVector(revisedBounds(1,j):revisedBounds(2,j)) = direction(revisedBounds(1,j):revisedBounds(2,j));
            uniqueDirs = cat(2,uniqueDirs,circ_mean(direction(revisedBounds(1,j):revisedBounds(2,j))'));
        end
        waveBounds = revisedBounds;
    end
end
end

%% additional functions
function [waveVector,boundries] = simpleWaveSegments(dataOrig,lengthLimit)
    waveVector = nan(1,length(dataOrig));
    data = dataOrig;
    data(isnan(data)) = 0;
    data(data~=0) = 1;
    data = find(data==0);
    dataEpochs = find(diff(data)>1);
    boundries = [data(dataEpochs)+1,data(dataEpochs+1)-1]';
    boundries(:,diff(boundries)<lengthLimit) = []; 
    for k = 1:size(boundries,2)
        waveVector(boundries(1,k):boundries(2,k)) = dataOrig(boundries(1,k):boundries(2,k));
    end
end

function [waveIndices] = getSegmentedWaves(data,boundaryIndices,devAngle)
data = reshape(data,[1,length(data)]);
dataFull = [nan data nan];
uniqueAngles = unique(dataFull);
uniqueAngles(isnan(uniqueAngles)) = [];

segWaveBounds = [];
inMat = [];
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