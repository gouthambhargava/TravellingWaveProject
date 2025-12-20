% load outputs and data
dataPath = 'F:\monkeyData\data';
gridType = 'Microelectrode';
subjectName='alpaH'; expDate = '210817'; protocolName = 'GRF_002'; 
sPos = 2; % spatial frequency: 0.5 (1), 1(2), 2 (3), 4 (4), 8 (5), all SFs (6). Note that the same code can be used for the size project also later where stimulus size is changed instead of spatial frequency
oriPos = 4; % orientation: 0 (1), 22.5 (2), 45 (3), 67.5 (4), 90 (5), 112.5 (6), 135 (7), 157.5 (8), all orientations (9)
[mData,goodElectrodes,timeVals,rfData,parameters] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sPos,oriPos);


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

waveVectorFull = cell(numFrequencyRanges,numTrials);
for i = 1:numTrials
    [waveVectorFull{1,i},~,~,~,waveType{1,i}] = getWaveParameters(outputs{1,i},squeeze(phases(:,i,:,1)),locList,timeVals,10,0, 4,boundryLims,2);
    [waveVectorFull{2,i},~,~,~,waveType{2,i}] = getWaveParameters(outputs{2,i},squeeze(phases(:,i,:,2)),locList,timeVals,10,0, 4,boundryLims,2);
end

intPts = nan(numTrials,length(timeVals));
overlap = 0.5;
for i = 1:numTrials
    [~,~,~,~,~,intPts(i,:)] = getOverlappingWaves(waveVector(i,:,1),waveBounds{1,i},waveVector(i,:,2),waveBounds{2,i},overlap);
end
boundries = [0.25 0.75];
%
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
colors = winter(3);

for j = 1:numTrials
    waveVectorSG = waveVectorFull{1,j};
    waveTypeSG = waveType{1,j};
    for i = 1:size(waveVectorSG,1)
        plot(timeVals,squeeze(waveVectorSG(i,:,:,1)).*trialId(j)-0.3,'LineWidth',1.2,'Color',colors(waveTypeSG(i),:))
        hold on
    end
    waveVectorFG = waveVectorFull{2,j};
    waveTypeFG = waveType{2,j};
    for i = 1:size(waveVectorFG,1)
        plot(timeVals,squeeze(waveVectorFG(i,:,:,1)).*trialId(j)-0.5,'LineWidth',1.2,'Color',colors(waveTypeFG(i),:))
        hold on
    end
end
yline(trialId,'--','LineWidth',0.5,'Color','black')
xlim(boundries)
ylim([0 numTrials+1])
yticks(0.5:1:(numTrials+0.5))
yticklabels(num2str(trialId'))
xlabel('Time (s)')
ylabel('Trials')
title('All classified waves:M1')

annotation('textbox',[0.91, 0.787, 0.1, 0.1],'String','Planar Wave','FontSize',12,'EdgeColor','none','Color',colors(1,:))
annotation('textbox',[0.91, 0.757, 0.1, 0.1],'String','Spiral Wave','FontSize',12,'EdgeColor','none','Color',colors(2,:))
annotation('textbox',[0.91, 0.727, 0.1, 0.1],'String','Complex Wave','FontSize',12,'EdgeColor','none','Color',colors(3,:))




