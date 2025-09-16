function [direction, pgd,spatialFreq] = getWaveMetrics(elecLocs,phaseData,burstLocs,waveMethod,neighbourLimit,numElectrodeCutoff)

% Inputs
% elecLocs is the list of electrode locations (x and y locations in the grid - dimensions reduced from 3D)
% phaseData is the instant phases of all elecs at the time point
% burstLocs is a logical vector, same length as the phaseData with gamma burst locations
% neighbourLimit - specifies the number of neighbours to be considered for
% calculating outputs for each electride
% waveMethod - set to 1 - wave metrics by considering all electrodes in the array.
%                     2 - wave metrics based on neighbouring clusters.
% Outputs described in circRegressEEG.
%% initialize the outputs
elecLocs = reshape(elecLocs,length(elecLocs),2);
numElectrodes = size(phaseData,1);
% if strcmpi(electrodeChoice,'all') %take all electrodes
%     burstElecs = 1:numGoodElectrodes;
% else
%     burstElecs = find(burstLocs);%find indices of electrodes that have bursts
% end

burstElecs = find(burstLocs);%find indices of electrodes that have bursts
 
if waveMethod ==1
    direction = nan;
    pgd = nan;
    spatialFreq = nan;
    if numel(burstElecs)>numElectrodeCutoff
        [direction,spatialFreq,~,~,pgd] = circRegMod(phaseData(burstElecs),elecLocs(burstElecs,:)); % change circRegressEEG to circRegMod
    end
        direction = repmat(direction,numElectrodes,1);
        spatialFreq = repmat(spatialFreq,numElectrodes,1);
        pgd = repmat(pgd,numElectrodes,1);
else
    direction = nan(numElectrodes,1);
    pgd = nan(numElectrodes,1);
    spatialFreq = nan(numElectrodes,1);
    distmat = squareform(pdist(elecLocs)); % define the distance matrix
    if numel(burstElecs)>numElectrodeCutoff
        % select electrode clusters based on distance
        for i = 1:length(elecLocs)
            [~,sortedIndices] = sort(distmat(i,:)); % sort the electrodes based on distance from current electrode in the loop
            neighbourElecs = sortedIndices(1:neighbourLimit);
            sigElecs = intersect(neighbourElecs,burstElecs);
            newLocs = elecLocs(sigElecs,:);
            newPhases = phaseData(sigElecs,:);
            if numel(newPhases)>2
                [direction(i),spatialFreq(i),~,~,pgd(i)] = circRegMod(newPhases,newLocs);
            end
        end
    end
end
end
