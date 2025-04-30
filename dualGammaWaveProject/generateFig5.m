%% generate fig 5
% load data for all ori combinations
% get bursts calculation for all
% generate outputs for both SG and FG for every ori combination
% run getWaveAndBurst overlap
% plot the results

% set up some parameters
minBurstSize = 100; % in ms
wobble = 5; % in deg
thresh = 0.5;
binEdges = 0:0.05:1;
electrodeFraction = 0.5;
electrodeChoice = 'selected';
arrayType = 'Microelectrode';
waveDetectionMethod = 1;

freqRangeList{1} = [25 35]; freqRangeList{2} = [40 50];
waveLengthLimit = 25;
sPos = 2; % spatial frequency: 0.5 (1), 1(2), 2 (3), 4 (4), 8 (5), all SFs (6). Note that the same code can be used for the size project also later where stimulus size is changed instead of spatial frequency
oriPos = 1:8; % orientation: 0 (1), 22.5 (2), 45 (3), 67.5 (4), 90 (5), 112.5 (6), 135 (7), 157.5 (8), all orientations (9)
stimPeriod = [0.25 0.75];
dataPath = 'G:\monkeyData\data';
gridType = 'Microelectrode';
%% load data
%for monkey 1
subjectName='alpaH'; expDate = '210817'; protocolName = 'GRF_002';

slowGammaOverlap1 = cell(1,length(oriPos));
fastGammaOverlap1 = cell(1,length(oriPos));
burstTS1 = cell(1,length(oriPos));
for i = 1:length(oriPos) % for all ori
    [allData,goodElectrodes,timeVals,rfData,parameters] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sPos,oriPos(i));

    %get electrode positions
    locList = zeros(length(goodElectrodes),2);
    gridLayout = rot90(reshape(1:81,[9,9]),2); %set the grid layout
    for gridi = 1:numel(goodElectrodes)
        [locList(gridi,1),locList(gridi,2)] = find(gridLayout==goodElectrodes(i));
    end

    burstMat = nan(size(allData));
    outputs = cell(length(freqRangeList),size(allData,2));
    for j = 1:size(allData,2) % for all trials
        for k = 1:length(freqRangeList)
            [burstMat(:,:,j),~,bandPhase] = getFilteredBurstsTW(squeeze(allData(:,j,:)),freqRangeList{k},[0.25 0.75],2,timeVals);
            outputs{k,j} = getTWCircParams(bandPhase,burstMat(:,:,j),timeVals,goodElectrodes,locList,electrodeFraction,electrodeChoice,arrayType,waveDetectionMethod);
        end
    end
    [slowGammaOverlap1{i},fastGammaOverlap1{i},burstTS1{i}] = getWaveAndBurstOverlap(burstTS,outputs,timeVals,burstLengthLimit,waveLengthLimit,waveWobble,binEdges);
end



%for monkey 2
subjectName='kesariH'; expDate = '270218'; protocolName = 'GRF_001';

slowGammaOverlap2 = cell(1,length(oriPos));
fastGammaOverlap2 = cell(1,length(oriPos));
burstTS2 = cell(1,length(oriPos));
for i = 1:length(oriPos) % for all ori
    [allData,goodElectrodes,timeVals,rfData,parameters] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sPos,oriPos(i));

    %get electrode positions
    locList = zeros(length(goodElectrodes),2);
    gridLayout = rot90(reshape(1:81,[9,9]),2); %set the grid layout
    for gridi = 1:numel(goodElectrodes)
        [locList(gridi,1),locList(gridi,2)] = find(gridLayout==goodElectrodes(i));
    end

    burstMat = nan(size(allData));
    outputs = cell(length(freqRangeList),size(allData,2));
    for j = 1:size(allData,2) % for all trials
        for k = 1:length(freqRangeList)
            [burstMat(:,:,j),~,bandPhase] = getFilteredBurstsTW(squeeze(allData(:,j,:)),freqRangeList{k},[0.25 0.75],2,timeVals);
            outputs{k,j} = getTWCircParams(bandPhase,burstMat(:,:,j),timeVals,goodElectrodes,locList,electrodeFraction,electrodeChoice,arrayType,waveDetectionMethod);
        end
    end
    [slowGammaOverlap2{i},fastGammaOverlap2{i},burstTS2{i}] = getWaveAndBurstOverlap(burstTS,outputs,timeVals,burstLengthLimit,waveLengthLimit,waveWobble,binEdges);
end



%% generate the plot
colorVals = cat(1,[52 148 186]./255,[236 112 22]./255);
subplot(4,4,[1,5,9])
burstData = sgBursts1{4};
burstData(burstData>0) = 1;
burstData(burstData==0) = nan;
multi = 1:size(burstData,1);
plot(binEdges(2:end),burstData.*multi','Color',colorVals(1,:))
% xlabel('Gamma Burst Bins')
ylabel(['Travelling Waves Locations:67.5',char(176)])
xlim([0.1 1])
ylim([1 size(burstData,1)])
title('Slow Gamma:M1')

subplot(4,4,[2,6,10])
burstData = fgBursts1{4};
burstData(burstData>0) = 1;
burstData(burstData==0) = nan;
multi = 1:size(burstData,1);
plot(binEdges(2:end),burstData.*multi','Color',colorVals(2,:))
% xlabel('Gamma Burst Bins')
ylabel(['Travelling Waves Locations:67.5',char(176)])
xlim([0.1 1])
ylim([1 size(burstData,1)])
title('Fast Gamma:M1')

subplot(4,4,[3,7,11])
burstData = sgBursts2{4};
burstData(burstData>0) = 1;
burstData(burstData==0) = nan;
multi = 1:size(burstData,1);
plot(binEdges(2:end),burstData.*multi','Color',colorVals(1,:))
% xlabel('Gamma Burst Bins')
ylabel(['Travelling Waves Locations:67.5',char(176)])
xlim([0.1 1])
ylim([1 size(burstData,1)])
title('Slow Gamma:M2')

subplot(4,4,[4,8,12])
burstData = fgBursts2{4};
burstData(burstData>0) = 1;
burstData(burstData==0) = nan;
multi = 1:size(burstData,1);
plot(binEdges(2:end),burstData.*multi','Color',colorVals(2,:))
% xlabel('Gamma Burst Bins')
ylabel(['Travelling Waves Locations:67.5',char(176)])
xlim([0.1 1])
ylim([1 size(burstData,1)])
title('Fast Gamma:M1')

subplot(4,4,13)
burstData = sum(sgBursts1{4});
plot(binEdges(2:end),burstData/max(burstData),'-o','LineWidth',1.2,'Color',colorVals(1,:))
hold on
burstData = sum(fgBursts1{4});
plot(binEdges(2:end),burstData/max(burstData),'-o','LineWidth',1.2,'Color',colorVals(2,:))
ylabel(['Ori-67.5',char(176)])
xlim([0.1 1])
title('TW distribution along \gamma bursts')

subplot(4,4,15)
burstData = sum(sgBursts2{4});
plot(binEdges(2:end),burstData/max(burstData),'-o','LineWidth',1.2,'Color',colorVals(1,:))
hold on
burstData = sum(fgBursts2{4});
plot(binEdges(2:end),burstData/max(burstData),'-o','LineWidth',1.2,'Color',colorVals(2,:))
xlim([0.1 1])
ylabel(['Ori-67.5',char(176)])
title('TW distribution along \gamma bursts')

subplot(4,4,14)
burstData = cell2mat(cellfun(@sum,sgBursts1','UniformOutput',false));
burstData = burstData./max(burstData,[],2);
burstDataErr = std(burstData);
plot(binEdges(2:end),mean(burstData),'-o','LineWidth',1.2,'Color',colorVals(1,:))
hold on
errorbar(binEdges(2:end),mean(burstData),burstDataErr,'LineWidth',1.2,'Color',colorVals(1,:))

burstData = cell2mat(cellfun(@sum,fgBursts1','UniformOutput',false));
burstData = burstData./max(burstData,[],2);
burstDataErr = std(burstData);
plot(binEdges(2:end),mean(burstData),'-o','LineWidth',1.2,'Color',colorVals(2,:))
hold on
errorbar(binEdges(2:end),mean(burstData),burstDataErr,'LineWidth',1.2,'Color',colorVals(2,:))
xlim([0.1 1])
xlabel('Gamma Burst Bins')
ylabel('All Ori')
title('TW distribution along \gamma bursts-all ori')

subplot(4,4,16)
burstData = cell2mat(cellfun(@sum,sgBursts2','UniformOutput',false));
burstData = burstData./max(burstData,[],2);
burstDataErr = std(burstData);
plot(binEdges(2:end),mean(burstData),'-o','LineWidth',1.2,'Color',colorVals(1,:))
hold on
errorbar(binEdges(2:end),mean(burstData),burstDataErr,'LineWidth',1.2,'Color',colorVals(1,:))

burstData = cell2mat(cellfun(@sum,fgBursts2','UniformOutput',false));
burstData = burstData./max(burstData,[],2);
burstDataErr = std(burstData);
plot(binEdges(2:end),mean(burstData),'-o','LineWidth',1.2,'Color',colorVals(2,:))
hold on
errorbar(binEdges(2:end),mean(burstData),burstDataErr,'LineWidth',1.2,'Color',colorVals(2,:))
xlim([0.1 1])
xlabel('Gamma Burst Bins')
ylabel('All Ori')
title('TW distribution along \gamma bursts-all ori')
%
annotation('textbox',[0.09,0.98, 0, 0], 'string', 'A','FontSize',20,'FontWeight','bold')
annotation('textbox',[0.3,0.98, 0, 0], 'string', 'B','FontSize',20,'FontWeight','bold')
annotation('textbox',[0.505,0.98, 0, 0], 'string', 'E','FontSize',20,'FontWeight','bold')
annotation('textbox',[0.72,0.98, 0, 0], 'string', 'F','FontSize',20,'FontWeight','bold')
annotation('textbox',[0.09,0.3, 0, 0], 'string', 'C','FontSize',20,'FontWeight','bold')
annotation('textbox',[0.3,0.3, 0, 0], 'string', 'D','FontSize',20,'FontWeight','bold')
annotation('textbox',[0.505,0.3, 0, 0], 'string', 'G','FontSize',20,'FontWeight','bold')
annotation('textbox',[0.72,0.3, 0, 0], 'string', 'H','FontSize',20,'FontWeight','bold')




