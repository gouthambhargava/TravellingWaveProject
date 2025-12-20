%% generate fig 5
% load data for all ori combinations
% get bursts calculation for all
% generate outputs for both SG and FG for every ori combination
% run getWaveAndBurst overlap
% plot the results

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
[X,Y] = meshgrid(1:1:9);
waveLengthLimit = 10;

thresholdFactor = 3;
% stimulusDurationS = [0 0.8]; % Stimulus duration to be highlighted
baselinePeriodS = [-0.5 0];
stimulusPeriodS = [0.25 0.75];
analysisPeriodS = [-0.5 1];
filterOrder = 4;

sPos = 2; % spatial frequency: 0.5 (1), 1(2), 2 (3), 4 (4), 8 (5), all SFs (6). Note that the same code can be used for the size project also later where stimulus size is changed instead of spatial frequency
oriPos = 1:8; % orientation: 0 (1), 22.5 (2), 45 (3), 67.5 (4), 90 (5), 112.5 (6), 135 (7), 157.5 (8), all orientations (9)
stimPeriod = [0.25 0.75];
dataPath = 'F:\monkeyData\data';
gridType = 'Microelectrode';

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
    [slowGammaOverlap1{i},fastGammaOverlap1{i},burstTS1{i}] = getWaveAndBurstOverlap(burstTS,outputs,timeVals,minBurstSize,waveLengthLimit,wobble,binEdges,goodElectrodes,thresh,3);
    [slowGammaOverlap1_5{i},fastGammaOverlap1_5{i},burstTS1{i}] = getWaveAndBurstOverlap(burstTS,outputs,timeVals,minBurstSize,waveLengthLimit,5,binEdges,goodElectrodes,thresh,2);
    [slowGammaOverlap1_10{i},fastGammaOverlap1_10{i},burstTS1{i}] = getWaveAndBurstOverlap(burstTS,outputs,timeVals,minBurstSize,waveLengthLimit,10,binEdges,goodElectrodes,thresh,2);

clear burstTS allData goodElectrodes filteredSignal
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
        [locList(gridi,1),locList(gridi,2)] = find(gridLayout==goodElectrodes(i));
    end

    % Do burst estimation
    for iFreq=1:numFrequencyRanges
        for iElec=1:numGoodElectrodes
            [~,~,~,burstTS(iElec,:,:,iFreq),filteredSignal(iElec,:,:,iFreq),~] = getHilbertBurst(squeeze(allData(iElec,:,:)),timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,freqRangeList{iFreq},filterOrder,1,analysisPeriodS);
        end
    end

    outputs = outputsTW2{i};
    [slowGammaOverlap2{i},fastGammaOverlap2{i},burstTS2{i}] = getWaveAndBurstOverlap(burstTS,outputs,timeVals,minBurstSize,waveLengthLimit,wobble,binEdges,goodElectrodes,thresh,3);
    [slowGammaOverlap2_5{i},fastGammaOverlap2_5{i},burstTS2{i}] = getWaveAndBurstOverlap(burstTS,outputs,timeVals,minBurstSize,waveLengthLimit,5,binEdges,goodElectrodes,thresh,2);
    [slowGammaOverlap2_10{i},fastGammaOverlap2_10{i},burstTS2{i}] = getWaveAndBurstOverlap(burstTS,outputs,timeVals,minBurstSize,waveLengthLimit,10,binEdges,goodElectrodes,thresh,2);
    
    clear burstTS allData goodElectrodes filteredSignal
end
toc


%% generate the plot
colorVals = cat(1,[52 148 186]./255,[236 112 22]./255);
subplot(5,4,[1,5])
burstData = slowGammaOverlap1{4};
burstData(burstData>0) = 1;
burstData(burstData==0) = nan;
multi = 1:size(burstData,1);
plot(binEdges(2:end),burstData.*multi','Color',colorVals(1,:))
% xlabel('Gamma Burst Bins')
ylabel(['Travelling Waves Locations:67.5',char(176)])
xlim([0.1 1])
ylim([1 size(burstData,1)])
title('Slow Gamma:M1')

subplot(5,4,[2,6])
burstData = fastGammaOverlap1{4};
burstData(burstData>0) = 1;
burstData(burstData==0) = nan;
multi = 1:size(burstData,1);
plot(binEdges(2:end),burstData.*multi','Color',colorVals(2,:))
% xlabel('Gamma Burst Bins')
ylabel(['Travelling Waves Locations:67.5',char(176)])
xlim([0.1 1])
ylim([1 size(burstData,1)])
title('Fast Gamma:M1')

subplot(5,4,[3,7])
burstData = slowGammaOverlap2{4};
burstData(burstData>0) = 1;
burstData(burstData==0) = nan;
multi = 1:size(burstData,1);
plot(binEdges(2:end),burstData.*multi','Color',colorVals(1,:))
% xlabel('Gamma Burst Bins')
ylabel(['Travelling Waves Locations:67.5',char(176)])
xlim([0.1 1])
ylim([1 size(burstData,1)])
title('Slow Gamma:M2')

subplot(5,4,[4,8])
burstData = fastGammaOverlap2{4};
burstData(burstData>0) = 1;
burstData(burstData==0) = nan;
multi = 1:size(burstData,1);
plot(binEdges(2:end),burstData.*multi','Color',colorVals(2,:))
% xlabel('Gamma Burst Bins')
ylabel(['Travelling Waves Locations:67.5',char(176)])
xlim([0.1 1])
ylim([1 size(burstData,1)])
title('Fast Gamma:M2')

subplot(5,4,9)
burstData = sum(slowGammaOverlap1{4});
plot(binEdges(2:end),burstData/max(burstData),'-o','LineWidth',1.2,'Color',colorVals(1,:))
hold on
burstData = sum(fastGammaOverlap1{4});
plot(binEdges(2:end),burstData/max(burstData),'-o','LineWidth',1.2,'Color',colorVals(2,:))
ylabel(['Ori-67.5',char(176)])
xlim([0.1 1])
title('TW distribution along \gamma bursts')

subplot(5,4,11)
burstData = sum(slowGammaOverlap2{4});
plot(binEdges(2:end),burstData/max(burstData),'-o','LineWidth',1.2,'Color',colorVals(1,:))
hold on
burstData = sum(fastGammaOverlap2{4});
plot(binEdges(2:end),burstData/max(burstData),'-o','LineWidth',1.2,'Color',colorVals(2,:))
xlim([0.1 1])
ylabel(['Ori-67.5',char(176)])
title('TW distribution along \gamma bursts')

subplot(5,4,10)
burstData = cell2mat(cellfun(@sum,slowGammaOverlap1','UniformOutput',false));
burstData = burstData./max(burstData,[],2);
burstDataErr = std(burstData);
plot(binEdges(2:end),mean(burstData),'-o','LineWidth',1.2,'Color',colorVals(1,:))
hold on
errorbar(binEdges(2:end),mean(burstData),burstDataErr,'LineWidth',1.2,'Color',colorVals(1,:))

burstData = cell2mat(cellfun(@sum,fastGammaOverlap1','UniformOutput',false));
burstData = burstData./max(burstData,[],2);
burstDataErr = std(burstData);
plot(binEdges(2:end),mean(burstData),'-o','LineWidth',1.2,'Color',colorVals(2,:))
hold on
errorbar(binEdges(2:end),mean(burstData),burstDataErr,'LineWidth',1.2,'Color',colorVals(2,:))
xlim([0.1 1])
xlabel('Gamma Burst Bins')
ylabel('All Ori')
title('TW distribution along \gamma bursts-all ori')

subplot(5,4,12)
burstData = cell2mat(cellfun(@sum,slowGammaOverlap2','UniformOutput',false));
burstData = burstData./max(burstData,[],2);
burstDataErr = std(burstData);
plot(binEdges(2:end),mean(burstData),'-o','LineWidth',1.2,'Color',colorVals(1,:))
hold on
errorbar(binEdges(2:end),mean(burstData),burstDataErr,'LineWidth',1.2,'Color',colorVals(1,:))

burstData = cell2mat(cellfun(@sum,fastGammaOverlap2','UniformOutput',false));
burstData = burstData./max(burstData,[],2);
burstDataErr = std(burstData);
plot(binEdges(2:end),mean(burstData),'-o','LineWidth',1.2,'Color',colorVals(2,:))
hold on
errorbar(binEdges(2:end),mean(burstData),burstDataErr,'LineWidth',1.2,'Color',colorVals(2,:))
xlim([0.1 1])
xlabel('Gamma Burst Bins')
ylabel('All Ori')
title('TW distribution along \gamma bursts-all ori')

% for 5 degree deviation
subplot(5,4,13)
burstData = sum(slowGammaOverlap1_5{4});
plot(binEdges(2:end),burstData/max(burstData),'-o','LineWidth',1.2,'Color',colorVals(1,:))
hold on
burstData = sum(fastGammaOverlap1_5{4});
plot(binEdges(2:end),burstData/max(burstData),'-o','LineWidth',1.2,'Color',colorVals(2,:))
ylabel(['Ori-67.5',char(176)])
xlim([0.1 1])
title('TW distribution along \gamma bursts')

subplot(5,4,15)
burstData = sum(slowGammaOverlap2_5{4});
plot(binEdges(2:end),burstData/max(burstData),'-o','LineWidth',1.2,'Color',colorVals(1,:))
hold on
burstData = sum(fastGammaOverlap2_5{4});
plot(binEdges(2:end),burstData/max(burstData),'-o','LineWidth',1.2,'Color',colorVals(2,:))
xlim([0.1 1])
ylabel(['Ori-67.5',char(176)])
title('TW distribution along \gamma bursts')

subplot(5,4,14)
burstData = cell2mat(cellfun(@sum,slowGammaOverlap1_5','UniformOutput',false));
burstData = burstData./max(burstData,[],2);
burstDataErr = std(burstData);
plot(binEdges(2:end),mean(burstData),'-o','LineWidth',1.2,'Color',colorVals(1,:))
hold on
errorbar(binEdges(2:end),mean(burstData),burstDataErr,'LineWidth',1.2,'Color',colorVals(1,:))

burstData = cell2mat(cellfun(@sum,fastGammaOverlap1_5','UniformOutput',false));
burstData = burstData./max(burstData,[],2);
burstDataErr = std(burstData);
plot(binEdges(2:end),mean(burstData),'-o','LineWidth',1.2,'Color',colorVals(2,:))
hold on
errorbar(binEdges(2:end),mean(burstData),burstDataErr,'LineWidth',1.2,'Color',colorVals(2,:))
xlim([0.1 1])
xlabel('Gamma Burst Bins')
ylabel('All Ori')
title('TW distribution along \gamma bursts-all ori')

subplot(5,4,16)
burstData = cell2mat(cellfun(@sum,slowGammaOverlap2_5','UniformOutput',false));
burstData = burstData./max(burstData,[],2);
burstDataErr = std(burstData);
plot(binEdges(2:end),mean(burstData),'-o','LineWidth',1.2,'Color',colorVals(1,:))
hold on
errorbar(binEdges(2:end),mean(burstData),burstDataErr,'LineWidth',1.2,'Color',colorVals(1,:))

burstData = cell2mat(cellfun(@sum,fastGammaOverlap2_5','UniformOutput',false));
burstData = burstData./max(burstData,[],2);
burstDataErr = std(burstData);
plot(binEdges(2:end),mean(burstData),'-o','LineWidth',1.2,'Color',colorVals(2,:))
hold on
errorbar(binEdges(2:end),mean(burstData),burstDataErr,'LineWidth',1.2,'Color',colorVals(2,:))
xlim([0.1 1])
xlabel('Gamma Burst Bins')
ylabel('All Ori')
title('TW distribution along \gamma bursts-all ori')

% for 10 degree deviation
subplot(5,4,17)
burstData = sum(slowGammaOverlap1_10{4});
plot(binEdges(2:end),burstData/max(burstData),'-o','LineWidth',1.2,'Color',colorVals(1,:))
hold on
burstData = sum(fastGammaOverlap1_10{4});
plot(binEdges(2:end),burstData/max(burstData),'-o','LineWidth',1.2,'Color',colorVals(2,:))
ylabel(['Ori-67.5',char(176)])
xlim([0.1 1])
title('TW distribution along \gamma bursts')

subplot(5,4,19)
burstData = sum(slowGammaOverlap2_10{4});
plot(binEdges(2:end),burstData/max(burstData),'-o','LineWidth',1.2,'Color',colorVals(1,:))
hold on
burstData = sum(fastGammaOverlap2_10{4});
plot(binEdges(2:end),burstData/max(burstData),'-o','LineWidth',1.2,'Color',colorVals(2,:))
xlim([0.1 1])
ylabel(['Ori-67.5',char(176)])
title('TW distribution along \gamma bursts')

subplot(5,4,18)
burstData = cell2mat(cellfun(@sum,slowGammaOverlap1_10','UniformOutput',false));
burstData = burstData./max(burstData,[],2);
burstDataErr = std(burstData);
plot(binEdges(2:end),mean(burstData),'-o','LineWidth',1.2,'Color',colorVals(1,:))
hold on
errorbar(binEdges(2:end),mean(burstData),burstDataErr,'LineWidth',1.2,'Color',colorVals(1,:))

burstData = cell2mat(cellfun(@sum,fastGammaOverlap1_10','UniformOutput',false));
burstData = burstData./max(burstData,[],2);
burstDataErr = std(burstData);
plot(binEdges(2:end),mean(burstData),'-o','LineWidth',1.2,'Color',colorVals(2,:))
hold on
errorbar(binEdges(2:end),mean(burstData),burstDataErr,'LineWidth',1.2,'Color',colorVals(2,:))
xlim([0.1 1])
xlabel('Gamma Burst Bins')
ylabel('All Ori')
title('TW distribution along \gamma bursts-all ori')

subplot(5,4,20)
burstData = cell2mat(cellfun(@sum,slowGammaOverlap2_10','UniformOutput',false));
burstData = burstData./max(burstData,[],2);
burstDataErr = std(burstData);
plot(binEdges(2:end),mean(burstData),'-o','LineWidth',1.2,'Color',colorVals(1,:))
hold on
errorbar(binEdges(2:end),mean(burstData),burstDataErr,'LineWidth',1.2,'Color',colorVals(1,:))

burstData = cell2mat(cellfun(@sum,fastGammaOverlap2_10','UniformOutput',false));
burstData = burstData./max(burstData,[],2);
burstDataErr = std(burstData);
plot(binEdges(2:end),mean(burstData),'-o','LineWidth',1.2,'Color',colorVals(2,:))
hold on
errorbar(binEdges(2:end),mean(burstData),burstDataErr,'LineWidth',1.2,'Color',colorVals(2,:))
xlim([0.1 1])
xlabel('Gamma Burst Bins')
ylabel('All Ori')
title('TW distribution along \gamma bursts-all ori')

annotation('textbox',[0.09,0.98, 0, 0], 'string', 'A','FontSize',20,'FontWeight','bold')
annotation('textbox',[0.3,0.98, 0, 0], 'string', 'B','FontSize',20,'FontWeight','bold')
annotation('textbox',[0.505,0.98, 0, 0], 'string', 'C','FontSize',20,'FontWeight','bold')
annotation('textbox',[0.72,0.98, 0, 0], 'string', 'D','FontSize',20,'FontWeight','bold')

annotation('textbox',[0.09,0.61, 0, 0], 'string', 'E','FontSize',20,'FontWeight','bold')
annotation('textbox',[0.3,0.61, 0, 0], 'string', 'F','FontSize',20,'FontWeight','bold')
annotation('textbox',[0.505,0.61, 0, 0], 'string', 'G','FontSize',20,'FontWeight','bold')
annotation('textbox',[0.72,0.61, 0, 0], 'string', 'H','FontSize',20,'FontWeight','bold')

annotation('textbox',[0.09,0.45, 0, 0], 'string', 'I','FontSize',20,'FontWeight','bold')
annotation('textbox',[0.3,0.45, 0, 0], 'string', 'J','FontSize',20,'FontWeight','bold')
annotation('textbox',[0.505,0.45, 0, 0], 'string', 'K','FontSize',20,'FontWeight','bold')
annotation('textbox',[0.72,0.45, 0, 0], 'string', 'L','FontSize',20,'FontWeight','bold')

annotation('textbox',[0.09,0.27, 0, 0], 'string', 'M','FontSize',20,'FontWeight','bold')
annotation('textbox',[0.3,0.27, 0, 0], 'string', 'N','FontSize',20,'FontWeight','bold')
annotation('textbox',[0.505,0.27, 0, 0], 'string', 'O','FontSize',20,'FontWeight','bold')
annotation('textbox',[0.72,0.27, 0, 0], 'string', 'P','FontSize',20,'FontWeight','bold')

annotation('textbox',[0.92,0.55, 0, 0], 'string', ['0',char(176)] ,'FontSize',12,'FontWeight','bold')
annotation('textbox',[0.92,0.38, 0, 0], 'string', ['5',char(176)],'FontSize',12,'FontWeight','bold')
annotation('textbox',[0.92,0.20, 0, 0], 'string', ['10',char(176)],'FontSize',12,'FontWeight','bold')

