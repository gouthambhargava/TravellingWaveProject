% load outputs and data
dataPath = 'F:\monkeyData\data';
gridType = 'Microelectrode';
subjectName='alpaH'; expDate = '210817'; protocolName = 'GRF_002'; 
sPos = 2; % spatial frequency: 0.5 (1), 1(2), 2 (3), 4 (4), 8 (5), all SFs (6). Note that the same code can be used for the size project also later where stimulus size is changed instead of spatial frequency
oriPos = 4; % orientation: 0 (1), 22.5 (2), 45 (3), 67.5 (4), 90 (5), 112.5 (6), 135 (7), 157.5 (8), all orientations (9)
[mData,goodElectrodes,timeVals] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sPos,oriPos);


load('alpaH_42_0.5T_selected_met2.mat')
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
        [~,~,~,burstTS(iElec,:,:,iFreq),filteredSignal(iElec,:,:,iFreq),~] = getHilbertBurst(squeeze(mData(iElec,:,:)),timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,freqRangeList{iFreq},filterOrder,1,analysisPeriodS);
        phases(iElec,:,:,iFreq) = angle(hilbert(squeeze(filteredSignal(iElec,:,:,iFreq))'))';    
    end
end

locList = zeros(length(goodElectrodes),2);
gridLayout = rot90(reshape(1:81,[9,9]),2); %set the grid layout
for gridi = 1:numGoodElectrodes
    [locList(gridi,1),locList(gridi,2)] = find(gridLayout==goodElectrodes(gridi));
end

waveVector1 = nan(numTrials,length(timeVals),numFrequencyRanges);
waveBounds1 = cell(numFrequencyRanges,numTrials);
directions = nan(numGoodElectrodes,length(timeVals),numTrials,numFrequencyRanges);

for i = 1:numFrequencyRanges
    for j = 1:numTrials
        [waveVector1(j,:,i),~,waveBounds1{i,j}] = getWaveSegments(outputs{i,j},timeVals,[],4,[0.25 0.75], 10);
        directions(:,:,j,i) = outputs{i,j}.direction;
    end
end

waveVectorFull1 = cell(numFrequencyRanges,numTrials);
for i = 1:numTrials
    [waveVectorFull1{1,i},~,~,~,waveType1{1,i}] = getWaveParameters(outputs{1,i},squeeze(phases(:,i,:,1)),locList,timeVals,10,0, 4,boundryLims,2);
    [waveVectorFull1{2,i},~,~,~,waveType1{2,i}] = getWaveParameters(outputs{2,i},squeeze(phases(:,i,:,2)),locList,timeVals,10,0, 4,boundryLims,2);
end

intPts1 = nan(numTrials,length(timeVals));
overlap = 0.5;
for i = 1:numTrials
    [~,~,~,~,~,intPts1(i,:)] = getOverlappingWaves(waveVector1(i,:,1),waveBounds1{1,i},waveVector1(i,:,2),waveBounds1{2,i},overlap);
end
numTrials1 = numTrials;

trialId1 = 1:numTrials1;

clear mData numTrials goodElectrodes rfData parameters outputs filteredSignal burstTS phases directions

%% for M2
% load outputs and data
dataPath = 'F:\monkeyData\data';
gridType = 'Microelectrode';
subjectName='kesariH'; expDate = '270218'; protocolName = 'GRF_001';
sPos = 2; % spatial frequency: 0.5 (1), 1(2), 2 (3), 4 (4), 8 (5), all SFs (6). Note that the same code can be used for the size project also later where stimulus size is changed instead of spatial frequency
oriPos = 4; % orientation: 0 (1), 22.5 (2), 45 (3), 67.5 (4), 90 (5), 112.5 (6), 135 (7), 157.5 (8), all orientations (9)
[mData,goodElectrodes,timeVals] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sPos,oriPos);


load('kesariH_42_0.5T_selected_met2.mat')
numTrials = size(outputs,2);

numGoodElectrodes = length(goodElectrodes);
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
        [~,~,~,burstTS(iElec,:,:,iFreq),filteredSignal(iElec,:,:,iFreq),~] = getHilbertBurst(squeeze(mData(iElec,:,:)),timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,freqRangeList{iFreq},filterOrder,1,analysisPeriodS);
        phases(iElec,:,:,iFreq) = angle(hilbert(squeeze(filteredSignal(iElec,:,:,iFreq))'))';    
    end
end

locList = zeros(length(goodElectrodes),2);
gridLayout = rot90(reshape(1:81,[9,9]),2); %set the grid layout
for gridi = 1:numGoodElectrodes
    [locList(gridi,1),locList(gridi,2)] = find(gridLayout==goodElectrodes(gridi));
end

waveVector2 = nan(numTrials,length(timeVals),numFrequencyRanges);
waveBounds2 = cell(numFrequencyRanges,numTrials);
directions = nan(numGoodElectrodes,length(timeVals),numTrials,numFrequencyRanges);

for i = 1:numFrequencyRanges
    for j = 1:numTrials
        [waveVector2(j,:,i),~,waveBounds2{i,j}] = getWaveSegments(outputs{i,j},timeVals,[],4,[0.25 0.75], 10);
        directions(:,:,j,i) = outputs{i,j}.direction;
    end
end

waveVectorFull2 = cell(numFrequencyRanges,numTrials);
for i = 1:numTrials
    [waveVectorFull2{1,i},~,~,~,waveType2{1,i}] = getWaveParameters(outputs{1,i},squeeze(phases(:,i,:,1)),locList,timeVals,10,0, 4,boundryLims,2);
    [waveVectorFull2{2,i},~,~,~,waveType2{2,i}] = getWaveParameters(outputs{2,i},squeeze(phases(:,i,:,2)),locList,timeVals,10,0, 4,boundryLims,2);
end

intPts2 = nan(numTrials,length(timeVals));
overlap = 0.5;
for i = 1:numTrials
    [~,~,~,~,~,intPts2(i,:)] = getOverlappingWaves(waveVector2(i,:,1),waveBounds2{1,i},waveVector2(i,:,2),waveBounds2{2,i},overlap);
end
numTrials2 = numTrials;
trialId2 = 1:numTrials2;

%%
boundries = [0.25 0.75];
subplot(1,4,1)
plotWavesAllTrials(waveVector1,timeVals,[0.25 0.75],intPts1)
ylim([0 36])
yticks(0.5:1:(numTrials1+0.5))
yticklabels(num2str(trialId1'))
xlabel('Time (s)')
ylabel('Trials')
title('All waves:M1')

subplot(1,4,2)
colors = winter(3);

for j = 1:numTrials1
    waveVectorSG = waveVectorFull1{1,j};
    waveTypeSG = waveType1{1,j};
    for i = 1:size(waveVectorSG,1)
        plot(timeVals,squeeze(waveVectorSG(i,:,:,1)).*trialId1(j)-0.3,'LineWidth',1.2,'Color',colors(waveTypeSG(i),:))
        hold on
    end
    waveVectorFG = waveVectorFull1{2,j};
    waveTypeFG = waveType1{2,j};
    for i = 1:size(waveVectorFG,1)
        plot(timeVals,squeeze(waveVectorFG(i,:,:,1)).*trialId1(j)-0.5,'LineWidth',1.2,'Color',colors(waveTypeFG(i),:))
        hold on
    end
end
yline(trialId2,'--','LineWidth',0.5,'Color','black')
xlim(boundries)
ylim([0 numTrials1+1])
yticks(0.5:1:(numTrials1+0.5))
yticklabels(num2str(trialId1'))
xlabel('Time (s)')
ylabel('Trials')
title('All classified waves:M1')

subplot(1,4,3)
plotWavesAllTrials(waveVector2,timeVals,[0.25 0.75],intPts2)
ylim([0 47])
yticks(0.5:1:(numTrials2+0.5))
yticklabels(num2str(trialId2'))
xlabel('Time (s)')
ylabel('Trials')
title('All waves:M2')

subplot(1,4,4)
colors = winter(3);

for j = 1:numTrials2
    waveVectorSG = waveVectorFull2{1,j};
    waveTypeSG = waveType2{1,j};
    for i = 1:size(waveVectorSG,1)
        plot(timeVals,squeeze(waveVectorSG(i,:,:,1)).*trialId2(j)-0.3,'LineWidth',1.2,'Color',colors(waveTypeSG(i),:))
        hold on
    end
    waveVectorFG = waveVectorFull2{2,j};
    waveTypeFG = waveType2{2,j};
    for i = 1:size(waveVectorFG,1)
        plot(timeVals,squeeze(waveVectorFG(i,:,:,1)).*trialId2(j)-0.5,'LineWidth',1.2,'Color',colors(waveTypeFG(i),:))
        hold on
    end
end
yline(trialId2,'--','LineWidth',0.5,'Color','black')
xlim(boundries)
ylim([0 numTrials2+1])
yticks(0.5:1:(numTrials2+0.5))
yticklabels(num2str(trialId2'))
xlabel('Time (s)')
ylabel('Trials')
title('All classified waves:M2')

annotation('textbox',[0.91, 0.787, 0.1, 0.1],'String','Planar Wave','FontSize',12,'EdgeColor','none','Color',colors(1,:))
annotation('textbox',[0.91, 0.757, 0.1, 0.1],'String','Spiral Wave','FontSize',12,'EdgeColor','none','Color',colors(2,:))
annotation('textbox',[0.91, 0.727, 0.1, 0.1],'String','Complex Wave','FontSize',12,'EdgeColor','none','Color',colors(3,:))

annotation('textbox',[0.095, 0.95, 0.02, 0.02], 'string', 'A','FontSize',24,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.305, 0.965, 0.02, 0.02], 'string', 'B','FontSize',24,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.51, 0.965, 0.02, 0.02], 'string', 'C','FontSize',24,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.72, 0.965, 0.02, 0.02], 'string', 'D','FontSize',24,'FontWeight','bold','EdgeColor','none')



