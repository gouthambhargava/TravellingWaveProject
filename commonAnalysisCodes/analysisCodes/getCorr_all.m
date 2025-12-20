%% get circular correlation (monte carlo method) for different methods
%load data
load('TravellingWaveProject\dualGammaWaveProject\data\alpaHM1.mat') %for alpaH
% load('TravellingWaveProject\dualGammaWaveProject\data\kesariHM1.mat') %for kesariH

% load the following for no threshold outputs
% load('TravellingWaveProject\dualGammaWaveProject\data\alpaHM1_noThresh.mat') %for alpaH
% load('TravellingWaveProject\dualGammaWaveProject\data\kesariHM1_noThresh.mat') %for kesariH

%% 1. with 5 degree mean variation across the full wave (used in paper)
% done for waves detected with gamma bursts and threshold fraction set at 50 for a 5 deg mean variance across waves 
% CC is done for mean wave angles
load('TravellingWaveProject\dualGammaWaveProject\data\alpaHM1.mat')
numFreqLimits = 2;
numTrials = size(outputs,2);
numTimePoints = length(timeVals);
waveVector = nan(numTimePoints,numTrials,numFreqLimits);
waveBounds = cell(numFreqLimits,numTrials);
newBounds = cell(1,numTrials);
allDirSg = nan(numTrials,numTimePoints);
allDirFg = nan(numTrials,numTimePoints);
allUniqueDirs = cell(1,numTrials);
emptyCell = nan(1,numTrials);

for i = 1:numFreqLimits
    for j = 1:numTrials
        [waveVector(:,j,i),uniqueDirs{i,j},waveBounds{i,j}] = getWaveSegments(outputs{i,j},timeVals,10,2,[0.25 0.75], 25);
    end
end

% get overlapping waves
for j = 1:numTrials
    [newBounds{j},allDirSg(j,:),allDirFg(j,:),allUniqueDirs{j},emptyCell(j)] = getOverlappingWaves(waveVector(:,j,1),waveBounds{1,j},waveVector(:,j,2),waveBounds{2,j},0.5);
end
numWaves1 = numel(cat(1,cell2mat(uniqueDirs(1,:)'),cell2mat(uniqueDirs(2,:)')))
allUniqueDirs(emptyCell==1) = [];
alpha = cell2mat(allUniqueDirs);
numel(alpha(1,:))
% get monte Carlo CC
[M1_corrCC,M1_corrP] = circCorrPermute(alpha(1,:), alpha(2,:), 999,1);

load('TravellingWaveProject\dualGammaWaveProject\data\kesariHM1.mat')
numFreqLimits = 2;
numTrials = size(outputs,2);
numTimePoints = length(timeVals);
waveVector = nan(numTimePoints,numTrials,numFreqLimits);
waveBounds = cell(numFreqLimits,numTrials);
newBounds = cell(1,numTrials);
allDirSg = nan(numTrials,numTimePoints);
allDirFg = nan(numTrials,numTimePoints);
allUniqueDirs = cell(1,numTrials);
emptyCell = nan(1,numTrials);

for i = 1:numFreqLimits
    for j = 1:numTrials
        [waveVector(:,j,i),uniqueDirs{i,j},waveBounds{i,j}] = getWaveSegments(outputs{i,j},timeVals,10,2,[0.25 0.75], 25);
    end
end
numWaves2 = numel(cat(1,cell2mat(uniqueDirs(1,:)'),cell2mat(uniqueDirs(2,:)')))

% get overlapping waves
for j = 1:numTrials
    [newBounds{j},allDirSg(j,:),allDirFg(j,:),allUniqueDirs{j},emptyCell(j)] = getOverlappingWaves(waveVector(:,j,1),waveBounds{1,j},waveVector(:,j,2),waveBounds{2,j},0.5);
end

allUniqueDirs(emptyCell==1) = [];
alpha = cell2mat(allUniqueDirs);
numel(alpha(1,:))
% get monte Carlo CC
[M2_corrCC,M2_corrP] = circCorrPermute(alpha(1,:), alpha(2,:), 999,1);

%% 2. with 5 degree variation across the whole wave
% done for waves detected with gamma bursts and threshold fraction set at
% 50 for a 5 deg variance across the entire wave
% CC is done for mean wave angles
load('TravellingWaveProject\dualGammaWaveProject\data\alpaHM1.mat')
numFreqLimits = 2;
numTrials = size(outputs,2);
numTimePoints = length(timeVals);
waveVector = nan(numTimePoints,numTrials,numFreqLimits);
waveBounds = cell(numFreqLimits,numTrials);
newBounds = cell(1,numTrials);
allDirSg = nan(numTrials,numTimePoints);
allDirFg = nan(numTrials,numTimePoints);
allUniqueDirs = cell(1,numTrials);
emptyCell = nan(1,numTrials);

for i = 1:numFreqLimits
    for j = 1:numTrials
        [waveVector(:,j,i),uniqueDirs{i,j},waveBounds{i,j}] = getWaveSegments(outputs{i,j},timeVals,5,3,[0.25 0.75], 25);
    end
end

% get overlapping waves
for j = 1:numTrials
    [newBounds{j},allDirSg(j,:),allDirFg(j,:),allUniqueDirs{j},emptyCell(j)] = getOverlappingWaves(waveVector(:,j,1),waveBounds{1,j},waveVector(:,j,2),waveBounds{2,j},0.5);
end

allUniqueDirs(emptyCell==1) = [];
alpha = cell2mat(allUniqueDirs);
size(alpha,2)
% get monte Carlo CC
[M1_corrCC,M1_corrP] = circCorrPermute(alpha(1,:), alpha(2,:), 999,1);

load('TravellingWaveProject\dualGammaWaveProject\data\kesariHM1.mat')
numFreqLimits = 2;
numTrials = size(outputs,2);
numTimePoints = length(timeVals);
waveVector = nan(numTimePoints,numTrials,numFreqLimits);
waveBounds = cell(numFreqLimits,numTrials);
newBounds = cell(1,numTrials);
allDirSg = nan(numTrials,numTimePoints);
allDirFg = nan(numTrials,numTimePoints);
allUniqueDirs = cell(1,numTrials);
emptyCell = nan(1,numTrials);

for i = 1:numFreqLimits
    for j = 1:numTrials
        [waveVector(:,j,i),uniqueDirs{i,j},waveBounds{i,j}] = getWaveSegments(outputs{i,j},timeVals,5,3,[0.25 0.75], 25);
    end
end

% get overlapping waves
for j = 1:numTrials
    [newBounds{j},allDirSg(j,:),allDirFg(j,:),allUniqueDirs{j},emptyCell(j)] = getOverlappingWaves(waveVector(:,j,1),waveBounds{1,j},waveVector(:,j,2),waveBounds{2,j},0.5);
end

allUniqueDirs(emptyCell==1) = [];
alpha = cell2mat(allUniqueDirs);
size(alpha,2)
% get monte Carlo CC
[M2_corrCC,M2_corrP] = circCorrPermute(alpha(1,:), alpha(2,:), 999,1);

%% 3. with 5 degree variation across the successive time points (to detect spiral waves) without gamma bursts or thresholding
% done for waves detected without gamma bursts and threshold fraction, for a 5 deg variance between each successive time point
% CC is done for all time points
load('TravellingWaveProject\dualGammaWaveProject\data\alpaHM1_noGamma.mat')
numFreqLimits = 2;
numTrials = size(outputs,2);
numTimePoints = length(timeVals);
waveVector = nan(numTimePoints,numTrials,numFreqLimits);
waveBounds = cell(numFreqLimits,numTrials);
newBounds = cell(1,numTrials);
allDirSg = nan(numTrials,numTimePoints);
allDirFg = nan(numTrials,numTimePoints);

for i = 1:numFreqLimits
    for j = 1:numTrials
        [waveVector(:,j,i),~,waveBounds{i,j}] = getWaveSegments(outputs{i,j},timeVals,5,2,[0.25 0.75], 25);
    end
end

% get overlapping waves
for j = 1:numTrials
    [newBounds{j},allDirSg(j,:),allDirFg(j,:),allUniqueDirs{j},emptyCell(j)] = getOverlappingWaves(waveVector(:,j,1),waveBounds{1,j},waveVector(:,j,2),waveBounds{2,j},0.5);
end

% get only overlapping time points
allDirSg(isnan(allDirFg)) = nan;
allDirFg(isnan(allDirSg)) = nan;
allDirSg(isnan(allDirSg)) = [];
allDirFg(isnan(allDirFg)) = [];

allUniqueDirs(emptyCell==1) = [];
alpha = cell2mat(allUniqueDirs);
% get monte Carlo CC
[corrCC,corrP] = circCorrPermute(alpha(1,:), alpha(2,:), 999,1)

load('TravellingWaveProject\dualGammaWaveProject\data\kesariHM1_noThresh.mat')
numFreqLimits = 2;
numTrials = size(outputs,2);
numTimePoints = length(timeVals);
waveVector = nan(numTimePoints,numTrials,numFreqLimits);
waveBounds = cell(numFreqLimits,numTrials);
newBounds = cell(1,numTrials);
allDirSg = nan(numTrials,numTimePoints);
allDirFg = nan(numTrials,numTimePoints);

for i = 1:numFreqLimits
    for j = 1:numTrials
        [waveVector(:,j,i),~,waveBounds{i,j}] = getWaveSegments(outputs{i,j},timeVals,5,2,[0.25 0.75], 25);
    end
end

% get overlapping waves
for j = 1:numTrials
    [newBounds{j},allDirSg(j,:),allDirFg(j,:),allUniqueDirs{j},emptyCell(j)] = getOverlappingWaves(waveVector(:,j,1),waveBounds{1,j},waveVector(:,j,2),waveBounds{2,j},0.5);
end

% get only overlapping time points
allDirSg(isnan(allDirFg)) = nan;
allDirFg(isnan(allDirSg)) = nan;
allDirSg(isnan(allDirSg)) = [];
allDirFg(isnan(allDirFg)) = [];

allUniqueDirs(emptyCell==1) = [];
alpha = cell2mat(allUniqueDirs);
% get monte Carlo CC
[corrCC,corrP] = circCorrPermute(alpha(1,:), alpha(2,:), 999,1)

%% 4.  with 5 degree variation across the successive time points (to detect spiral waves) without gamma bursts or thresholding
% done for waves detected without gamma bursts and threshold fraction, for a 5 deg variance between each successive time point
% waves have been classified into planar and spiral and CC done seperately for each
% CC is done for mean waves
% load data
sPos = 2;
oriPos = 4;
req = 2;
dataPath = 'F:\monkeyData\data';
gridType = 'Microelectrode';
freqRangeList{1} = [20 35]; freqRangeList{2} = [40 60];
subjectName='alpaH'; expDate = '210817'; protocolName = 'GRF_002'; 
[allData,goodElectrodes,timeVals,~,~] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sPos,oriPos);

% get grid layout 
numFrequencyRanges = 2;
numChans = 1:numel(goodElectrodes);
numGoodElectrodes = numel(numChans);
freqRangeList{1} = [20 35]; freqRangeList{2} = [40 60];
numTrials = size(outputs,2);
locList = zeros(length(numChans),2);
gridLayout = rot90(reshape(1:81,[9,9]),2); %set the grid layout
for gridi = 1:numel(numChans)
    [locList(gridi,1),locList(gridi,2)] = find(gridLayout==goodElectrodes(gridi));
end

% define outputs
numTimePoints = length(timeVals);
waveVector = nan(numTimePoints,numTrials,numFreqLimits);
waveBounds = cell(numFreqLimits,numTrials);
newBounds = cell(1,numTrials);
allDirSg = nan(numTrials,numTimePoints);
allDirFg = nan(numTrials,numTimePoints);
filteredSignal = zeros(numel(goodElectrodes),numTrials,numTimePoints,numFrequencyRanges);
directionSG = nan(1,numTimePoints,numTrials);
directionFG = nan(1,numTimePoints,numTrials);
phaseSG = zeros(numGoodElectrodes,numTimePoints,numFrequencyRanges);
phaseFG = zeros(numGoodElectrodes,numTimePoints,numFrequencyRanges);
waveTypeSG = cell(1,numTrials);
waveTypeFG = cell(1,numTrials);

% Do burst estimation
thresholdFactor = 3;
baselinePeriodS = [-0.5 0];
stimulusPeriodS = [0.25 0.75];
analysisPeriodS = [-0.5 1];
filterOrder = 4;
for iFreq=1:numFrequencyRanges
    for iElec=1:numGoodElectrodes
        [~,~,~,~,filteredSignal(iElec,:,:,iFreq),~] = getHilbertBurst(squeeze(allData(iElec,:,:)),timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,freqRangeList{iFreq},filterOrder,1,analysisPeriodS);
    end
end

for i = 1:numTrials
    directionSG(:,:,i) = outputs{1,i}.direction';
    directionFG(:,:,i) = outputs{2,i}.direction';
end

directionSG = repmat(directionSG,numGoodElectrodes,1);
directionFG = repmat(directionFG,numGoodElectrodes,1);

for i = 1:numTrials
    phaseSG(:,:,i) = angle(hilbert(squeeze(filteredSignal(:,i,:,1))'))';
    phaseFG(:,:,i) = angle(hilbert(squeeze(filteredSignal(:,i,:,2))'))';
end

% calculate the wave segments
for i = 1:numFrequencyRanges
    for j = 1:numTrials
        [waveVector(:,j,i),~,waveBounds{i,j}] = getWaveSegments(outputs{i,j},timeVals,5,2,[0.25 0.75], 25);
    end
end

% get wave and phase frames
for i = 1:numTrials
    [~,~,waveTypeSG{i}] = getMultiWaveType(phaseSG(:,:,i),directionSG(:,:,i),locList,waveBounds{1,i},2);
    [~,~,waveTypeFG{i}] = getMultiWaveType(phaseFG(:,:,i),directionFG(:,:,i),locList,waveBounds{2,i},2);
end

% seperate out the planar and the spiral waves
planarBoundsSG = cell(1,numTrials);
spiralBoundsSG = cell(1,numTrials);
planarBoundsFG = cell(1,numTrials);
spiralBoundsFG = cell(1,numTrials);

for i = 1:numTrials
    planarWaveId = find(waveTypeSG{i}==1);
    planarBoundsSG{i} = waveBounds{1,i}(:,planarWaveId);
    spiralWaveId = find(waveTypeSG{i}>1);
    spiralBoundsSG{i} = waveBounds{1,i}(:,spiralWaveId);
end

for i = 1:numTrials
    planarWaveId = find(waveTypeFG{i}==1);
    planarBoundsFG{i} = waveBounds{2,i}(:,planarWaveId);
    spiralWaveId = find(waveTypeFG{i}>1);
    spiralBoundsFG{i} = waveBounds{2,i}(:,spiralWaveId);
end

% find overlapping waves 
%all overlapping waves
for j = 1:numTrials
    [~,allDirSg(j,:),allDirFg(j,:)] = getOverlappingWaves(waveVector(:,j,1),waveBounds{1,j},waveVector(:,j,2),waveBounds{2,j},0.5);
end

%planar overlapping waves
newBoundsPlanar = cell(1,numTrials);
allUniqueDirsP = cell(1,numTrials);
for j = 1:numTrials
    [newBoundsPlanar{j},~,~,allUniqueDirsP{j},~] = getOverlappingWaves(waveVector(:,j,1),planarBoundsSG{j},waveVector(:,j,2),planarBoundsFG{j},0.5);
end

%spiral overlapping waves
newBoundsSpiral = cell(1,numTrials);
allUniqueDirsS = cell(1,numTrials);
for j = 1:numTrials
    [newBoundsSpiral{j},~,~,allUniqueDirsS{j},~] = getOverlappingWaves(waveVector(:,j,1),spiralBoundsSG{j},waveVector(:,j,2),spiralBoundsFG{j},0.5);
end

%CC for averaged wavepoint
alpha1 = cell2mat(allUniqueDirsP);
alpha1(:,isnan(alpha1(1,:))) = []; 
alpha2 = cell2mat(allUniqueDirsS);
alpha2(:,isnan(alpha2(1,:))) = []; 

[corrCCPlanar,corrPPlanar] = circCorrPermute(alpha1(1,:), alpha1(2,:), 999,1);
[corrCCSpiral,corrPSpiral] = circCorrPermute(alpha2(1,:), alpha2(2,:), 999,1);

