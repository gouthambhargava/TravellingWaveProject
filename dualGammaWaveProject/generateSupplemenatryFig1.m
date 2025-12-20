% generate supplementary figure 2
% set up some parameters
% load waves for all ori 
% outputsTW1 = {outputsTWA21,outputsTWA22,outputsTWA23,outputsTWA24,outputsTWA25,outputsTWA26,outputsTWA27,outputsTWA28};
% outputsTW2 = {outputsTWK21,outputsTWK22,outputsTWK23,outputsTWK24,outputsTWK25,outputsTWK26,outputsTWK27,outputsTWK28};

% parameters for wave detection
wobble1 = 5; % in deg
wobble2 = 10;
lengthLimit = 10;
freqRangeList{1} = [20 35]; freqRangeList{2} = [40 60];

numFrequencyRanges = numel(freqRangeList);

%% get wave data for 5 and 10 deg for default ori and SF
% for M1
load('alpaH_42_0.5T_selected_met1.mat')
numTrials = length(outputs);
numTrials1 = numTrials;
waveVector1_5 = nan(numTrials,length(timeVals),numFrequencyRanges);
waveBounds1_5 = cell(numFrequencyRanges,numTrials);
waveVector1_10 = nan(numTrials,length(timeVals),numFrequencyRanges);
waveBounds1_10 = cell(numFrequencyRanges,numTrials);
dirSG1_5 = nan(numTrials,length(timeVals));
dirFG1_5 = nan(numTrials,length(timeVals));
dirSG1_10 = nan(numTrials,length(timeVals));
dirFG1_10 = nan(numTrials,length(timeVals));

for iFreq=1:numFrequencyRanges
        for iTrials=1:numTrials
            [waveVector1_5(iTrials,:,iFreq),~,waveBounds1_5{iFreq,iTrials}] = getWaveSegments(outputs{iFreq,iTrials},timeVals,wobble1,2,[0.25 0.75], lengthLimit);
            [waveVector1_10(iTrials,:,iFreq),~,waveBounds1_10{iFreq,iTrials}] = getWaveSegments(outputs{iFreq,iTrials},timeVals,wobble2,2,[0.25 0.75], lengthLimit);
        end
end

for i = 1:numTrials
    [~,dirSG1_5(i,:),dirFG1_5(i,:)] = getOverlappingWaves(waveVector1_5(i,:,1),waveBounds1_5{1,i},waveVector1_5(i,:,2),waveBounds1_5{2,i},0.5);
    [~,dirSG1_10(i,:),dirFG1_10(i,:)] = getOverlappingWaves(waveVector1_10(i,:,1),waveBounds1_10{1,i},waveVector1_10(i,:,2),waveBounds1_10{2,i},0.5);
end

ovWaves1_5 = dirSG1_5;
ovWaves1_5(isnan(dirFG1_5)) = nan;
ovWaves1_5(~isnan(ovWaves1_5)) = 1;
ovWaves1_10 = dirSG1_10;
ovWaves1_10(isnan(dirFG1_10)) = nan;
ovWaves1_10(~isnan(ovWaves1_10)) = 1;
clear dirSG1_5 dirSG1_10 dirFG1_5 dirFG1_10 waveBounds1_5 waveBounds1_10 i numTrials iFreq iTrails outputs

% for M2
load('kesariH_42_0.5T_selected_met1.mat')
numTrials = length(outputs);
numTrials2 = numTrials;
waveVector2_5 = nan(numTrials,length(timeVals),numFrequencyRanges);
waveBounds2_5 = cell(numFrequencyRanges,numTrials);
waveVector2_10 = nan(numTrials,length(timeVals),numFrequencyRanges);
waveBounds2_10 = cell(numFrequencyRanges,numTrials);
dirSG2_5 = nan(numTrials,length(timeVals));
dirFG2_5 = nan(numTrials,length(timeVals));
dirSG2_10 = nan(numTrials,length(timeVals));
dirFG2_10 = nan(numTrials,length(timeVals));

for iFreq=1:numFrequencyRanges
        for iTrials=1:numTrials
            [waveVector2_5(iTrials,:,iFreq),~,waveBounds2_5{iFreq,iTrials}] = getWaveSegments(outputs{iFreq,iTrials},timeVals,wobble1,2,[0.25 0.75], lengthLimit);
            [waveVector2_10(iTrials,:,iFreq),~,waveBounds2_10{iFreq,iTrials}] = getWaveSegments(outputs{iFreq,iTrials},timeVals,wobble2,2,[0.25 0.75], lengthLimit);
        end
end

for i = 1:numTrials
    [~,dirSG2_5(i,:),dirFG2_5(i,:)] = getOverlappingWaves(waveVector2_5(i,:,1),waveBounds2_5{1,i},waveVector2_5(i,:,2),waveBounds2_5{2,i},0.5);
    [~,dirSG2_10(i,:),dirFG2_10(i,:)] = getOverlappingWaves(waveVector2_10(i,:,1),waveBounds2_10{1,i},waveVector2_10(i,:,2),waveBounds2_10{2,i},0.5);
end

ovWaves2_5 = dirSG2_5;
ovWaves2_5(isnan(dirFG2_5)) = nan;
ovWaves2_5(~isnan(ovWaves2_5)) = 1;
ovWaves2_10 = dirSG2_10;
ovWaves2_10(isnan(dirFG2_10)) = nan;
ovWaves2_10(~isnan(ovWaves2_10)) = 1;
clear dirSG2_5 dirSG2_10 dirFG2_5 dirFG2_10 waveBounds2_5 waveBounds2_10 i numTrials iFreq iTrails


%% plot the data
colorVals = cat(1,[52 148 186]./255,[236 112 22]./255);
subplot(3,4,[1,5,9])
trialsNum1 = 0.6:numTrials1;
plotWavesAllTrials(waveVector1_5,timeVals,[0.25 0.75])
hold on
plot(timeVals,ovWaves1_5.*trialsNum1','LineWidth',0.8,'Color','black')
ylim([0 numTrials1+0.5])
ylabel('Trials')
xlabel('Time (s)')
title(['All waves- 5',char(176),' deviation:M1'])

subplot(3,4,[2,6,10])
plotWavesAllTrials(waveVector1_10,timeVals,[0.25 0.75])
hold on
plot(timeVals,ovWaves1_10.*trialsNum1','LineWidth',0.8,'Color','black')
ylim([0 numTrials1+0.5])
ylabel('Trials')
xlabel('Time (s)')
title(['All waves- 10',char(176),' deviation:M1'])

subplot(3,4,[3,7,11])
trialsNum2 = 0.6:numTrials2;
plotWavesAllTrials(waveVector2_5,timeVals,[0.25 0.75])
hold on
plot(timeVals,ovWaves2_5.*trialsNum2','LineWidth',0.8,'Color','black')
ylim([0 numTrials2+0.5])
ylabel('Trials')
xlabel('Time (s)')
title(['All waves- 5',char(176),' deviation:M2'])


subplot(3,4,[4,8,12])
plotWavesAllTrials(waveVector2_10,timeVals,[0.25 0.75])
hold on
plot(timeVals,ovWaves2_10.*trialsNum2','LineWidth',0.8,'Color','black')
ylim([0 numTrials2+0.5])
ylabel('Trials')
xlabel('Time (s)')
title(['All waves- 10',char(176),' deviation:M2'])
annotation('textbox',[0.09,0.98, 0, 0], 'string', 'A','FontSize',20,'FontWeight','bold')
annotation('textbox',[0.3,0.98, 0, 0], 'string', 'B','FontSize',20,'FontWeight','bold')
annotation('textbox',[0.505,0.98, 0, 0], 'string', 'C','FontSize',20,'FontWeight','bold')
annotation('textbox',[0.72,0.98, 0, 0], 'string', 'D','FontSize',20,'FontWeight','bold')

annotation('textbox',[0.01, 0.587, 0.1, 0.1],'String','Slow Gamma','FontSize',12,'Color',colorVals(1,:),'EdgeColor','none')
annotation('textbox',[0.01, 0.550, 0.1, 0.1],'String','Fast Gamma','FontSize',12,'Color',colorVals(2,:),'EdgeColor','none')
annotation('textbox',[0.01, 0.510, 0.1, 0.1],'String','Overlapping Waves','FontSize',12,'Color','black','EdgeColor','none')
