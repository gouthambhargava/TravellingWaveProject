function [direction, pgd,spatialFreq] = getWaveMetrics(locList,phiVals,burstLocs,waveMethod,neighbourLimit)

% Inputs
% locList is the list of electrode locations (x and y locations in the grid - dimensions reduced from 3D)
% phiVals is the instant phases of all elecs at the time point
% burstLocs is a logical vector, same length as the phiVals with gamma burst locations
% neighbourLimit - specifies the number of neighbours to be considered for
% calculating outputs for each electride
% waveMethod - set to 1 - wave metrics by considering all electrodes in the array.
%                     2 - wave metrics based on neighbouring clusters.
% Outputs described in circRegressEEG.
%% initialize the outputs
locList = reshape(locList,length(locList),2);
numElectrodes = size(phiVals,1);
burstElecs = find(burstLocs);%find indices of electrodes that have bursts

if waveMethod ==1

    [direction,spatialFreq,~,~,pgd] = circRegMod(phiVals(burstElecs),locList(burstElecs,:)); % change circRegressEEG to circRegMod
    direction = repmat(wrapTo2Pi(direction),numElectrodes,1);
    spatialFreq = repmat(spatialFreq,numElectrodes,1);
    pgd = repmat(pgd,numElectrodes,1);
else
    direction = nan(numElectrodes,1);
    pgd = nan(numElectrodes,1);
    spatialFreq = nan(numElectrodes,1);
    distmat = squareform(pdist(locList)); % define the distance matrix
    % select electrode clusters based on distance
    for i = 1:length(locList)
        [~,sortedIndices] = sort(distmat(i,:)); % sort the electrodes based on distance from current electrode in the loop
        neighbourElecs = sortedIndices(1:neighbourLimit);
        sigElecs = intersect(neighbourElecs,burstElecs);
        newLocs = locList(sigElecs,:);
        newPhases = phiVals(sigElecs,:);
        if numel(newPhases)>2
            [direction(i),spatialFreq(i),~,~,pgd(i)] = circRegMod(newPhases,newLocs);
        end
    end
end
end
