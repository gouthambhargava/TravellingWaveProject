function [slowGammaOverlap,fastGammaOverlap,burstTS] = getBurstOverlap(outputsTW,sPos,oriPos,thresh,minBurstSize,binEdges,req)
% function only for generating figure 5. Identifies overlapping gamma
% bursts and subsequently the variation of wave locations within these
% overlapping bursts
%% load data
dataPath = 'F:\monkeyData\data';
gridType = 'Microelectrode';
% sPos = 2; % spatial frequency: 0.5 (1), 1(2), 2 (3), 4 (4), 8 (5), all SFs (6). Note that the same code can be used for the size project also later where stimulus size is changed instead of spatial frequency
% oriPos = 4;% orientation: 0 (1), 22.5 (2), 45 (3), 67.5 (4), 90 (5), 112.5 (6), 135 (7), 157.5 (8), all orientations (9)

if req ==1
    subjectName='alpaH'; expDate = '210817'; protocolName = 'GRF_002'; 
    [mData,goodElectrodes,timeVals,~,~] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sPos,oriPos);
else
    % load data for M2
    subjectName='kesariH'; expDate = '270218'; protocolName = 'GRF_001';  
    [mData,goodElectrodes,timeVals,~,~] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sPos,oriPos);
end
% burst calculation parameters
thresholdFactor = 3;
baselinePeriodS = [-0.5 0];
stimulusPeriodS = [0.25 0.75];
analysisPeriodS = [0.25 0.75];
filterOrder = 4;
freqRangeList{1} = [20 35]; freqRangeList{2} = [40 60];
numFrequencyRanges = numel(freqRangeList);

for iFreq=1:numFrequencyRanges
    for iElec=1:numel(goodElectrodes)
        [~,~,~,burstTS(iElec,:,:,iFreq),~] = getHilbertBurst(squeeze(mData(iElec,:,:)),timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,freqRangeList{iFreq},filterOrder,1,analysisPeriodS);
     end
end
clear analysisPeriodS req baselinePeriodS dataPath expDate filterOrder freqRangeList gridType iElec iFreq oriPos sPos stimulusPeriodS subjectName thresholdFactor 
%% generate the binned data
% thresh = 0.5;
burstTS(isnan(burstTS)) = 0;
elecFrac = numel(goodElectrodes)*thresh;
burstFrac = squeeze(sum(burstTS));
burstFrac(burstFrac<elecFrac) = 0;
burstFrac(burstFrac==0) = nan;
% minBurstSize = 100;
lims = [0.25 0.75];
% wobble = 5;
[allBoundsSG,~] = getSegmentedBursts(burstFrac(:,:,1),lims,timeVals,minBurstSize);
[allBoundsFG,~] = getSegmentedBurst(burstFrac(:,:,2),lims,timeVals,minBurstSize);

allDir1 = getWaveSegments(outputsTW(1,:),timeVals,2);
allDir2 = getWaveSegments(outputsTW(2,:),timeVals,2);

allDir = cat(3,allDir1,allDir2);
allDir(~isnan(allDir)) = 1;
allDir(isnan(allDir)) = 0;
clear allDir1 allDir2

slowGammaOverlap = [];
fastGammaOverlap = [];
% binEdges = 0:0.1:1;

% get dir values corrosponding to the bursts
for i = 1:numel(allBoundsSG)
    for j = 1:size(allBoundsSG{i},2)
        int = allDir(i,allBoundsSG{i}(1,j):allBoundsSG{i}(2,j),1);
        int1 = allBoundsSG{i}(1,j):allBoundsSG{i}(2,j);
        int1 = (int1-min(int1))/(max(int1)-min(int1));
        int1(int==0) = [];
        y = histcounts(int1,binEdges);
        slowGammaOverlap = cat(1,slowGammaOverlap,y);
    end
    for j = 1:size(allBoundsFG{i},2)
        int = allDir(i,allBoundsFG{i}(1,j):allBoundsFG{i}(2,j),2);
        int1 = allBoundsFG{i}(1,j):allBoundsFG{i}(2,j);
        int1 = (int1-min(int1))/(max(int1)-min(int1));
        int1(int==0) = [];
        y = histcounts(int1,binEdges);
        fastGammaOverlap = cat(1,fastGammaOverlap,y);
    end
end
clear int int1 i j y
end

function [allBounds,allBursts] = getSegmentedBurst(bursts,Lims,timeVals,minBurstSize)
% set boundryLims as a vector of times outside of which direction vectors
% will not be considered
% if timeFlag = 1, append time values, else append the indices instead
boundryLims = [dsearchn(timeVals',Lims(1)),dsearchn(timeVals',Lims(2))];
bursts(:,1:boundryLims(1)-1) = nan;
bursts(:,boundryLims(2)+1:length(timeVals)) = nan;
%% find the travelling wave indices for all trials and plot
allBounds = cell(1,size(bursts,1));
allBursts = nan(size(bursts));
    for i = 1:size(bursts,1)
        boundries = [];
        burstInt = zeros(1,length(timeVals));
        burstInt(~isnan(bursts(i,:))) = 1;
        burstData = find(burstInt==0);
        burstEpochs = find(diff(burstData)>1);
        boundries = cat(2,boundries,[burstData(burstEpochs)+1;burstData(burstEpochs+1)-1]); 
        boundries(:,diff(boundries)<5) = [];
        boundries(:,diff(boundries)<minBurstSize) = [];
        for k = 1:size(boundries,2)
            allBursts(i,boundries(1,k):boundries(2,k)) = 1;
        end
        allBounds{i} = boundries;
    end
end
