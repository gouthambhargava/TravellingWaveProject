% get mean, SD of direction distributions
% outputs = outputsTWA24; % change the output condition to get different results
numTrials = size(outputs,2);
numFrequencyRanges = size(outputs,1);
timePts = size(outputs{1,1}.direction,2);
waveVector = nan(numTrials,timePts,numFrequencyRanges);
waveBounds = cell(numFrequencyRanges,numTrials);
uniqueDirs = cell(numFrequencyRanges,numTrials);
wobble = 0;
segOption = 3;
lengthLimit = 10;
overlap = 0.5;
boundryLims = [0.25 0.75];
for i = 1:numTrials
    for j = 1:numFrequencyRanges
        [waveVector(i,:,j),uniqueDirs{j,i},waveBounds{j,i}] = getWaveSegments(outputs{j,i},timeVals,wobble,segOption,boundryLims, lengthLimit);
    end
end

allDirSg = nan(numTrials,timePts);
allDirFg = nan(numTrials,timePts);
uniqueDirsOl = cell(1,numTrials);
for k = 1:numTrials
    [~,allDirSg(k,:),allDirFg(k,:),uniqueDirsOl{k}] = getOverlappingWaves(waveVector(k,:,1),waveBounds{1,k},waveVector(k,:,2),waveBounds{2,k},overlap);
end

% mean, median and SD of nonoverlapping directions
sgAnglesTemp = waveVector(:,:,1);
sgAnglesTemp(isnan(sgAnglesTemp)) = [];
fgAnglesTemp = waveVector(:,:,2);
fgAnglesTemp(isnan(fgAnglesTemp)) = [];
direcStatsSG = wrapTo360(rad2deg([circ_mean(sgAnglesTemp),circ_median(sgAnglesTemp),circ_std(sgAnglesTemp)]));
direcStatsFG = wrapTo360(rad2deg([circ_mean(fgAnglesTemp),circ_median(fgAnglesTemp),circ_std(fgAnglesTemp)]));

% mean, median and SD of overlapping directions
allDirSg(isnan(allDirSg)) = [];
allDirFg(isnan(allDirFg)) = [];
direcStatsOvSG = wrapTo360(rad2deg([circ_mean(allDirSg),circ_median(allDirSg),circ_std(allDirSg)]));
direcStatsOvFG = wrapTo360(rad2deg([circ_mean(allDirFg),circ_median(allDirFg),circ_std(allDirFg)]));

% get CC and pvals
[~,rhoVal,pVal,rhoVal2,pVal2] = simpleOvWrapper(outputs,timeVals,wobble,segOption,lengthLimit,overlap, 1,1);

%% get speed and duration (of waves and gamma bursts) and compare
% outputs = outputsTW1{4};

numTrials = size(outputs,2);
numFrequencyRanges = size(outputs,1);
timePts = size(outputs{1,1}.direction,2);
waveBounds = cell(numFrequencyRanges,numTrials);
burstBounds = cell(numFrequencyRanges,numTrials);
wobble = 0;
segOption = 3;
lengthLimit = 10;
boundryLims = [0.25 0.75];

freqRangeList{1} = [20 35]; freqRangeList{2} = [40 60];
numGoodElectrodes = length(goodElectrodes);
numFrequencyRanges = numel(freqRangeList);

% Do burst estimation
thresholdFactor = 3;
baselinePeriodS = [-0.5 0];
stimulusPeriodS = [0.25 0.75];
analysisPeriodS = [-0.5 1];
filterOrder = 4;
for iFreq=1:numFrequencyRanges
    burstLengthFull = [];
    for iElec=1:numGoodElectrodes
        burstLength = getHilbertBurst(squeeze(allData(iElec,:,:)),timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,freqRangeList{iFreq},filterOrder,1,analysisPeriodS);
        burstLengthFull = cat(2,burstLengthFull,cell2mat(burstLength));
    end
    allBursts{iFreq} = burstLengthFull;
end


for i = 1:numTrials
    for j = 1:numFrequencyRanges
        [~,~,waveBounds{j,i}] = getWaveSegments(outputs{j,i},timeVals,wobble,segOption,boundryLims, lengthLimit);    
        speed(i,:,j) = outputs{j,i}.speed(1,:);

    end
end




allDurSG = diff(cell2mat(waveBounds(1,:)),[],1)*0.5;
allDurFG = diff(cell2mat(waveBounds(2,:)),[],1)*0.5;
meanDur = [mean(allDurSG),mean(allDurFG)]; % in milliseconds
medianDur = [median(allDurSG),median(allDurFG)]; % in milliseconds

stdDur = [std(allDurSG),std(allDurFG)];
[p, ~, stats] = ranksum(allDurSG,allDurFG);


allSpeedSG = speed(:,:,1);
allSpeedSG(isnan(allSpeedSG)) = [];
allSpeedFG = speed(:,:,2);
allSpeedFG(isnan(allSpeedFG)) = [];
meanSpeed = [mean(allSpeedSG),mean(allSpeedFG)]; % in m/s
stdSpeed = [std(allSpeedSG),std(allSpeedFG)];

%% get overlap and CC for method 2
%load alpa and kesari HM1 and HM2
load('D:\IISC_work\TWGitScripts\TravellingWaveProject\dualGammaWaveProject\data\kesariHM2.mat')
numTrials = size(outputs,2);
numFrequencyRanges = size(outputs,1);
segOption = 4;
wobble = 0;
lengthLimit = 15;
boundryLims = [0.25 0.75];
overlap = 0.25;
% get wave segments
for i = 1:numTrials
    for j = 1:numFrequencyRanges
        [waveVector(i,:,j),~,waveBounds{j,i}] = getWaveSegments(outputs{j,i},timeVals,wobble,segOption,boundryLims, lengthLimit);    
        allDir(:,:,i,j) = outputs{j,i}.direction;
    end
end
% get overlapping points
for k = 1:numTrials
    [ovBounds{k},~,~,~,emptyTrials(k),intPts(k,:)] = getOverlappingWaves(waveVector(k,:,1),waveBounds{1,k},waveVector(k,:,2),waveBounds{2,k},overlap);
end
% get original directions for overlapping points
ovBounds(emptyTrials==1) = [];
allDir(:,:,emptyTrials==1,:) = [];

% method 1 - mean of full wave
allOverlaps = [];
for i = 1:numel(find(emptyTrials==0))
    bounds = ovBounds{i};
    for j = 1:size(bounds{1},2)
        overlaps(:,j,1) = circMeanNan(allDir(:,bounds{1}(1,j):bounds{1}(2,j),1),2);
        overlaps(:,j,2) = circMeanNan(allDir(:,bounds{2}(1,j):bounds{2}(2,j),2),2);
    end
    allOverlaps = cat(2,allOverlaps,overlaps);
end

% method 2 - mean of only overlapping segments
% intPts(emptyTrials==1,:) = [];
% allOverlaps = [];
% for i = 1:numel(find(emptyTrials==0))
%     [~,bounds] = simpleWaveSegments(intPts(i,:),0);
%     for j = 1:size(bounds,2)
%         overlaps(:,j,1) = circMeanNan(allDir(:,bounds(1,j):bounds(2,j),1),2);
%         overlaps(:,j,2) = circMeanNan(allDir(:,bounds(1,j):bounds(2,j),2),2);
%     end
%     allOverlaps = cat(2,allOverlaps,overlaps);
% end

allOverlaps(allOverlaps==0) = nan;
dir1 = allOverlaps(:,:,1);
dir2 = allOverlaps(:,:,2);
dir1(isnan(dir2)) = nan;
dir2(isnan(dir1)) = nan;
dir1(isnan(dir1)) = [];
dir2(isnan(dir2)) = [];


[cc(1), cc(2)] = circCorrPermute(dir1,dir2,50,1);
% [cc(1), cc(2)] = circ_corrcc(dir1,dir2);
clear allDir numTrials emptyTrials intPts allOverlaps overlaps waveVector waveBounds bounds
%% get kruskal wallis/anova for gamma bins
dataPath = 'F:\monkeyData\data';
gridType = 'Microelectrode';
% set up some parameters
minBurstSize = 25; % in ms
wobble = 0; % in degww
thresh = 0.5;
binEdges = 0:0.05:1;
electrodeFraction = 0.5;
electrodeChoice = 'selected';
arrayType = 'Microelectrode';
waveDetectionMethod = 1;

freqRangeList{1} = [20 35]; freqRangeList{2} = [40 60];
numFrequencyRanges = numel(freqRangeList);
% [X,Y] = meshgrid(1:1:9);
waveLengthLimit = 10;

thresholdFactor = 3;
% stimulusDurationS = [0 0.8]; % Stimulus duration to be highlighted
baselinePeriodS = [-0.5 0];
stimulusPeriodS = [0.25 0.75];
analysisPeriodS = [-0.5 1];
filterOrder = 4;

sPos = 2; % spatial frequency: 0.5 (1), 1(2), 2 (3), 4 (4), 8 (5), all SFs (6). Note that the same code can be used for the size project also later where stimulus size is changed instead of spatial frequency
oriPos = 1:8; % orientation: 0 (1), 22.5 (2), 45 (3), 67.5 (4), 90 (5), 112.5 (6), 135 (7), 157.5 (8), all orientations (9)
% stimPeriod = [0.25 0.75];


%for monkey 1
subjectName='alpaH'; expDate = '210817'; protocolName = 'GRF_002';

slowGammaOverlap1 = cell(1,length(oriPos));
fastGammaOverlap1 = cell(1,length(oriPos));
slowGammaOverlap1_5 = cell(1,length(oriPos));
fastGammaOverlap1_5 = cell(1,length(oriPos));
slowGammaOverlap1_10 = cell(1,length(oriPos));
fastGammaOverlap1_10 = cell(1,length(oriPos));
burstTS1 = cell(1,length(oriPos));
tic
for i = 1:length(oriPos) % for all ori
    [allData,goodElectrodes,timeVals] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sPos,oriPos(i));
    numGoodElectrodes = length(goodElectrodes);
    numTrials = size(allData,2);
    burstTS = nan(numGoodElectrodes,numTrials,length(timeVals),numFrequencyRanges);
    %get electrode positions
    locList = zeros(length(goodElectrodes),2);
    gridLayout = rot90(reshape(1:81,[9,9]),2); %set the grid layout
    for gridi = 1:numel(goodElectrodes)
        [locList(gridi,1),locList(gridi,2)] = find(gridLayout==goodElectrodes(gridi));
    end


    % Do burst estimation
    for iFreq=1:numFrequencyRanges
        for iElec=1:numGoodElectrodes
            [~,~,~,burstTS(iElec,:,:,iFreq),~,~] = getHilbertBurst(squeeze(allData(iElec,:,:)),timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,freqRangeList{iFreq},filterOrder,1,analysisPeriodS);
        end
    end
    outputs = outputsTW1{i};
    [slowGammaOverlap1{i},fastGammaOverlap1{i},~] = getWaveAndBurstOverlap(burstTS,outputs,timeVals,minBurstSize,waveLengthLimit,wobble,binEdges,goodElectrodes,thresh,3);
    [slowGammaOverlap1_5{i},fastGammaOverlap1_5{i},~] = getWaveAndBurstOverlap(burstTS,outputs,timeVals,minBurstSize,waveLengthLimit,5,binEdges,goodElectrodes,thresh,2);
    [slowGammaOverlap1_10{i},fastGammaOverlap1_10{i},~] = getWaveAndBurstOverlap(burstTS,outputs,timeVals,minBurstSize,waveLengthLimit,10,binEdges,goodElectrodes,thresh,2);

clear burstTS allData goodElectrodes filteredSignal outputs 
end
toc

% for monkey 2
tic
subjectName='kesariH'; expDate = '270218'; protocolName = 'GRF_001';

slowGammaOverlap2 = cell(1,length(oriPos));
fastGammaOverlap2 = cell(1,length(oriPos));
slowGammaOverlap2_5 = cell(1,length(oriPos));
fastGammaOverlap2_5 = cell(1,length(oriPos));
slowGammaOverlap2_10 = cell(1,length(oriPos));
fastGammaOverlap2_10 = cell(1,length(oriPos));
burstTS2 = cell(1,length(oriPos));
for i = 1:length(oriPos) % for all ori
    [allData,goodElectrodes,timeVals,rfData,parameters] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sPos,oriPos(i));
    numGoodElectrodes = length(goodElectrodes);
    numTrials = size(allData,2);
    burstTS = nan(numGoodElectrodes,numTrials,length(timeVals),numFrequencyRanges);

    %get electrode positions
    locList = zeros(length(goodElectrodes),2);
    gridLayout = rot90(reshape(1:81,[9,9]),2); %set the grid layout
    for gridi = 1:numel(goodElectrodes)
        [locList(gridi,1),locList(gridi,2)] = find(gridLayout==goodElectrodes(gridi));
    end

    % Do burst estimation
    for iFreq=1:numFrequencyRanges
        for iElec=1:numGoodElectrodes
            [~,~,~,burstTS(iElec,:,:,iFreq),filteredSignal(iElec,:,:,iFreq),~] = getHilbertBurst(squeeze(allData(iElec,:,:)),timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,freqRangeList{iFreq},filterOrder,1,analysisPeriodS);
        end
    end

    outputs = outputsTW2{i};
    [slowGammaOverlap2{i},fastGammaOverlap2{i},~] = getWaveAndBurstOverlap(burstTS,outputs,timeVals,minBurstSize,waveLengthLimit,wobble,binEdges,goodElectrodes,thresh,3);
    [slowGammaOverlap2_5{i},fastGammaOverlap2_5{i},~] = getWaveAndBurstOverlap(burstTS,outputs,timeVals,minBurstSize,waveLengthLimit,5,binEdges,goodElectrodes,thresh,2);
    [slowGammaOverlap2_10{i},fastGammaOverlap2_10{i},~] = getWaveAndBurstOverlap(burstTS,outputs,timeVals,minBurstSize,waveLengthLimit,10,binEdges,goodElectrodes,thresh,2);
    
    clear burstTS allData goodElectrodes filteredSignal outputs
end
toc


%% get stats
%% for M1
slowGammaOverlapAll = {slowGammaOverlap1,slowGammaOverlap1_5,slowGammaOverlap1_10};
fastGammaOverlapAll = {fastGammaOverlap1,fastGammaOverlap1_5,fastGammaOverlap1_10};

for j = 1:numel(slowGammaOverlapAll)
    group = cell2mat(slowGammaOverlapAll{j}');
    group(sum(group,2)==0,:) = [];
    numSamples = size(group,2);
    multiplier = 1:numSamples;
    idx = [];
    allBins = [];
    for i = 1:numSamples
        idx = cat(2,idx,multiplier(i)*ones(1,size(group,1)));
        allBins = cat(1,allBins,group(:,i));
    end
    [SGM1p(j),tbl1{j},stats] = kruskalwallis(allBins,idx,"off");
    SGM1c(:,:,j) = multcompare(stats);
end

for j = 1:numel(fastGammaOverlapAll)
    group = cell2mat(fastGammaOverlapAll{j}');
    group(sum(group,2)==0,:) = [];
    numSamples = size(group,2);
    multiplier = 1:numSamples;
    idx = [];
    allBins = [];
    for i = 1:numSamples
        idx = cat(2,idx,multiplier(i)*ones(1,size(group,1)));
        allBins = cat(1,allBins,group(:,i));
    end
    [FGM1p(j),tbl2{j},stats] = kruskalwallis(allBins,idx,"off");
    FGM1c(:,:,j) = multcompare(stats);
end

%% for M2
slowGammaOverlapAll = {slowGammaOverlap2,slowGammaOverlap2_5,slowGammaOverlap2_10};
fastGammaOverlapAll = {fastGammaOverlap2,fastGammaOverlap2_5,fastGammaOverlap2_10};

for j = 1:numel(slowGammaOverlapAll)
    group = cell2mat(slowGammaOverlapAll{j}');
    group(sum(group,2)==0,:) = [];
    numSamples = size(group,2);
    multiplier = 1:numSamples;
    idx = [];
    allBins = [];
    for i = 1:numSamples
        idx = cat(2,idx,multiplier(i)*ones(1,size(group,1)));
        allBins = cat(1,allBins,group(:,i));
    end
    [SGM2p(j),tbl1_1{j},stats] = kruskalwallis(allBins,idx,"off");
    SGM2c(:,:,j) = multcompare(stats);
end

for j = 1:numel(fastGammaOverlapAll)
    group = cell2mat(fastGammaOverlapAll{j}');
    group(sum(group,2)==0,:) = [];
    numSamples = size(group,2);
    multiplier = 1:numSamples;
    idx = [];
    allBins = [];
    for i = 1:numSamples
        idx = cat(2,idx,multiplier(i)*ones(1,size(group,1)));
        allBins = cat(1,allBins,group(:,i));
    end
    [FGM2p(j),tbl2_2{j},stats2] = kruskalwallis(allBins,idx,"off");
    FGM2c(:,:,j) = multcompare(stats2);
end
% clear group i j numSamples idx allBins stats slowGammaOverlapAll slowGammaOverlap1 slowGammaOverlap1_5 slowGammaOverlap1_10 fastGammaOverlapAll fastGammaOverlap1 fastGammaOverlap1_5 fastGammaOverlap1_10
% clear slowGammaOverlap2 slowGammaOverlap2_5 slowGammaOverlap2_10 fastGammaOverlap2 fastGammaOverlap2_5 fastGammaOverlap2_10
% clear dataPath gridType minBurstSize wobble thres binEdges electrodeFraction electrodeChoice arrayType waveDetectionMethod freqRangeList numFrequencyRanges X Y waveLengthLimit thresholdFactor baselinePeriodS stimulusPeriodS analysisPeriodS filterOrder sPos oriPos stimPeriod