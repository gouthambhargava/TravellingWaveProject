% load outputs and data
load('D:\IISC_work\gitScripts\FinalWobbleScripts\alpa24Data.mat')
load('D:\IISC_work\TWGitScripts\TravellingWaveProject\dualGammaWaveProject\data\alpaHM2.mat')

numTrials = size(outputs,2);
wobble = 0;
segOption = 4;
lengthLimit = 10;
boundryLims = [0.25 0.75];

freqRangeList{1} = [20 35]; freqRangeList{2} = [40 60];
numGoodElectrodes = length(goodElectrodes);
numFrequencyRanges = numel(freqRangeList);
burstTS = nan(numGoodElectrodes,numTrials,length(timeVals),numFrequencyRanges);
filteredSignal = nan(numGoodElectrodes,numTrials,length(timeVals),numFrequencyRanges);
phases = nan(numGoodElectrodes,numTrials,length(timeVals),numFrequencyRanges);

% Do burst estimation
thresholdFactor = 3;
baselinePeriodS = [-0.5 0];
stimulusPeriodS = [0.25 0.75];
analysisPeriodS = [-0.5 1];
filterOrder = 4;
for iFreq=1:numFrequencyRanges
    for iElec=1:numGoodElectrodes
        [~,~,~,burstTS(iElec,:,:,iFreq),filteredSignal(iElec,:,:,iFreq),~] = getHilbertBurst(squeeze(allData(iElec,:,:)),timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,freqRangeList{iFreq},filterOrder,1,analysisPeriodS);
        phases(iElec,:,:,iFreq) = angle(hilbert(squeeze(filteredSignal(iElec,:,:,iFreq))'))';    
    end
end

locList = zeros(length(goodElectrodes),2);
gridLayout = rot90(reshape(1:81,[9,9]),2); %set the grid layout
for gridi = 1:numGoodElectrodes
    [locList(gridi,1),locList(gridi,2)] = find(gridLayout==goodElectrodes(gridi));
end

waveVector = nan(numTrials,length(timeVals),numFrequencyRanges);
waveBounds = cell(numFrequencyRanges,numTrials);
directions = nan(numGoodElectrodes,length(timeVals),numTrials,numFrequencyRanges);

for i = 1:numFrequencyRanges
    for j = 1:numTrials
        [waveVector(j,:,i),~,waveBounds{i,j}] = getWaveSegments(outputs{i,j},timeVals,[],4,[0.25 0.75], 10);
        directions(:,:,j,i) = outputs{i,j}.direction;
    end
end

waveTypeSG = cell(1,numTrials);
waveTypeFG = cell(1,numTrials);

for i = 1:numTrials
    [~, ~,waveTypeSG{i},~,~,~] = interpAndWaveClassi(squeeze(phases(:,i,:,1)),directions(:,:,i,1),locList,waveBounds{1,i},1,1);
    [~, ~,waveTypeFG{i},~,~,~] = interpAndWaveClassi(squeeze(phases(:,i,:,2)),directions(:,:,i,2),locList,waveBounds{2,i},1,1);
end

waveVectorFull = nan(3,numTrials,length(timeVals),numFrequencyRanges);

for i = 1:numTrials
    waveType = waveTypeSG{i};
    bounds = waveBounds{1,i};
    for j = 1:size(bounds,2)
        waveVectorFull(waveType(j),i,bounds(1,j):bounds(2,j),1) = 1;
    end
end

for i = 1:numTrials
    waveType = waveTypeFG{i};
    bounds = waveBounds{2,i};
    for j = 1:size(bounds,2)
        waveVectorFull(waveType(j),i,bounds(1,j):bounds(2,j),2) = 1;
    end
end


intPts = nan(numTrials,length(timeVals));
overlap = 0.5;
for i = 1:numTrials
    [~,~,~,~,~,intPts(i,:)] = getOverlappingWaves(waveVector(i,:,1),waveBounds{1,i},waveVector(i,:,2),waveBounds{2,i},overlap);
end
boundries = [0.25 0.75];

trialId = 1:numTrials;
subplot(1,2,1)
plotWavesAllTrials(waveVector,timeVals,[0.25 0.75],intPts)
ylim([0 36])
yticks(0.5:1:(numTrials+0.5))
yticklabels(num2str(trialId'))
xlabel('Time (s)')
ylabel('Trials')
title('All waves:M1')

subplot(1,2,2)
colors = parula(3);
for i = 1:size(waveVectorFull,1)
    plot(timeVals,squeeze(waveVectorFull(i,:,:,1)).*trialId'-0.3,'LineWidth',1.2,'Color',colors(i,:))
    hold on
end
for i = 1:size(waveVectorFull,1)
    plot(timeVals,squeeze(waveVectorFull(i,:,:,2)).*trialId'-0.5,'LineWidth',1.2,'Color',colors(i,:))
    hold on
end
yline(trialId,'--','LineWidth',0.5,'Color','black')
xlim(boundries)
ylim([0 36])
yticks(0.5:1:(numTrials+0.5))
yticklabels(num2str(trialId'))
xlabel('Time (s)')
ylabel('Trials')
title('All classified waves:M1')


plot(timeVals,overlap.*trialId'-0.4,'LineWidth',1,'Color','black')


