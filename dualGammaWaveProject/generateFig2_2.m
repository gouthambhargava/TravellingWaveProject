%% generate fig 2
% load data
dataPath = 'F:\monkeyData\data';
gridType = 'Microelectrode';
subjectName='alpaH'; expDate = '210817'; protocolName = 'GRF_002'; 
sPos = 2; % spatial frequency: 0.5 (1), 1(2), 2 (3), 4 (4), 8 (5), all SFs (6). Note that the same code can be used for the size project also later where stimulus size is changed instead of spatial frequency
oriPos = 4; % orientation: 0 (1), 22.5 (2), 45 (3), 67.5 (4), 90 (5), 112.5 (6), 135 (7), 157.5 (8), all orientations (9)
[~,goodElectrodes1,timeVals] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sPos,oriPos);
load('alpaH_42_0.5T_selected_met1.mat');
outputsTW1 = outputs;

load('kesariH_42_0.5T_selected_met1.mat');
subjectName='kesariH'; expDate = '270218'; protocolName = 'GRF_001';
[~,goodElectrodes2] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sPos,oriPos);
outputsTW2 = outputs;

%set params
numTrials1 = size(outputsTW1,2);
numTrials2 = size(outputsTW2,2);

%% define some parameters for wave detection
lengthLimit = 10; %ms
boundryLims = [0.25 0.75];
wobbleLim = 0; %degree
segOption = 3;
numFrequencyRanges = size(outputsTW1,1);

% segment waves and get overlapping waves for all trials
%initialize outputs
waveVector1 = nan(numTrials1,length(timeVals),numFrequencyRanges);
waveBounds1 = cell(numFrequencyRanges,numTrials1);

for i = 1:numTrials1
    for j = 1:numFrequencyRanges
        [waveVector1(i,:,j),~,waveBounds1{j,i}] = getWaveSegments(outputsTW1{j,i},timeVals,wobbleLim,segOption,boundryLims, lengthLimit);
    end
end

% find overlapping waves 
intPts1 = nan(numTrials1,length(timeVals));
overlap = 0.5;
for i = 1:numTrials1
    [~,~,~,~,~,intPts1(i,:)] = getOverlappingWaves(waveVector1(i,:,1),waveBounds1{1,i},waveVector1(i,:,2),waveBounds1{2,i},overlap);
end

% get wave data for M2
waveVector2 = nan(numTrials2,length(timeVals),numFrequencyRanges);
waveBounds2 = cell(numFrequencyRanges,numTrials2);

for i = 1:numTrials2
    for j = 1:numFrequencyRanges
        [waveVector2(i,:,j),~,waveBounds2{j,i}] = getWaveSegments(outputsTW2{j,i},timeVals,wobbleLim,segOption,boundryLims, lengthLimit);
    end
end

% find overlapping waves 
intPts2 = nan(numTrials2,length(timeVals));
overlap = 0.5;
for i = 1:numTrials2
    [~,~,~,~,~,intPts2(i,:)] = getOverlappingWaves(waveVector2(i,:,1),waveBounds2{1,i},waveVector2(i,:,2),waveBounds2{2,i},overlap);
end

% plot all trials with overlapping waves        
trialsNum1 = 0.6:numTrials1;
trialLabels1 = 1:numTrials1;

trialsNum2 = 0.6:numTrials2;
trialLabels2 = 1:numTrials2;

%% plot
subplot(1,2,1)
plotWavesAllTrials(waveVector1,timeVals,[0.25 0.75])
hold on
plot(timeVals,intPts1.*trialsNum1','LineWidth',0.8,'Color','black')
hold on
ylim([0 numTrials1])
ylabel('Trials')
xlabel('Time (s)')
yticks(0.5:1:numTrials1)
yticklabels(num2str(trialLabels1'))
title('Overlapping waves across all trials:M1')

subplot(1,2,2)
plotWavesAllTrials(waveVector2,timeVals,[0.25 0.75])
hold on
plot(timeVals,intPts2.*trialsNum2','LineWidth',0.8,'Color','black')
hold on
ylim([0 numTrials2])
ylabel('Trials')
xlabel('Time (s)')
yticks(0.5:1:numTrials2)
yticklabels(num2str(trialLabels2'))
title('Overlapping waves across all trials:M2')


colorVals = cat(1,[52 148 186]./255,[236 112  22]./255);

annotation('textbox',[0.05,0.95, 0, 0], 'string', 'A','FontSize',30,'FontWeight','bold')
annotation('textbox',[0.5,0.95, 0, 0], 'string', 'B','FontSize',30,'FontWeight','bold')
annotation('textbox',[0.01, 0.587, 0.1, 0.1],'String','Slow Gamma','FontSize',12,'Color',colorVals(1,:),'EdgeColor','none')
annotation('textbox',[0.01, 0.550, 0.1, 0.1],'String','Fast Gamma','FontSize',12,'Color',colorVals(2,:),'EdgeColor','none')
annotation('textbox',[0.01, 0.510, 0.1, 0.1],'String','Overlapping Waves','FontSize',12,'Color','black','EdgeColor','none')



