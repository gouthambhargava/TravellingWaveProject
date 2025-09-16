%generate fig 4 - final
% get data
% load data for M1
% generate outputs for every ori combination
% plot


% set up some parameters
wobble = 5; % in deg
thresh = 0.5;
electrodeFraction = 0.5;
electrodeChoice = 'selected';
arrayType = 'Microelectrode';
waveDetectionMethod = 1;
orientations = [0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5];
segmentationMethod = 3;

freqRangeList{1} = [25 35]; freqRangeList{2} = [40 50];
waveLengthLimit = 25;% in ms
sPos = 2; % spatial frequency: 0.5 (1), 1(2), 2 (3), 4 (4), 8 (5), all SFs (6). Note that the same code can be used for the size project also later where stimulus size is changed instead of spatial frequency
oriPos = 1:8; % orientation: 0 (1), 22.5 (2), 45 (3), 67.5 (4), 90 (5), 112.5 (6), 135 (7), 157.5 (8), all orientations (9)
stimPeriod = [0.25 0.75];
dataPath = 'G:\monkeyData\data';
gridType = 'Microelectrode';
%% load data and get wave parameters
%for monkey 1
subjectName='alpaH'; expDate = '210817'; protocolName = 'GRF_002';

for i = 1:length(oriPos) % for all ori
    [allData,goodElectrodes,timeVals,rfData,parameters] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sPos,oriPos(i));

    %get electrode positions
    locList = zeros(length(goodElectrodes),2);
    gridLayout = rot90(reshape(1:81,[9,9]),2); %set the grid layout
    for gridi = 1:numel(goodElectrodes)
        [locList(gridi,1),locList(gridi,2)] = find(gridLayout==goodElectrodes(i));
    end

    burstMat = nan(size(allData));
    outputsTW1 = cell(length(freqRangeList),size(allData,2));
    for j = 1:size(allData,2) % for all trials
        for k = 1:length(freqRangeList)
            [burstMat(:,:,j),~,bandPhase] = getFilteredBurstsTW(squeeze(allData(:,j,:)),freqRangeList{k},[0.25 0.75],2,timeVals);
            outputsTW1{k,j} = getTWCircParams(bandPhase,burstMat(:,:,j),timeVals,goodElectrodes,locList,electrodeFraction,electrodeChoice,arrayType,waveDetectionMethod);
        end
    end
end

%for monkey 2
subjectName='kesariH'; expDate = '270218'; protocolName = 'GRF_001';

for i = 1:length(oriPos) % for all ori
    [allData,goodElectrodes,timeVals,rfData,parameters] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sPos,oriPos(i));

    %get electrode positions
    locList = zeros(length(goodElectrodes),2);
    gridLayout = rot90(reshape(1:81,[9,9]),2); %set the grid layout
    for gridi = 1:numel(goodElectrodes)
        [locList(gridi,1),locList(gridi,2)] = find(gridLayout==goodElectrodes(i));
    end

    burstMat = nan(size(allData));
    outputsTW2 = cell(length(freqRangeList),size(allData,2));
    for j = 1:size(allData,2) % for all trials
        for k = 1:length(freqRangeList)
            [burstMat(:,:,j),~,bandPhase] = getFilteredBurstsTW(squeeze(allData(:,j,:)),freqRangeList{k},[0.25 0.75],2,timeVals);
            outputsTW2{k,j} = getTWCircParams(bandPhase,burstMat(:,:,j),timeVals,goodElectrodes,locList,electrodeFraction,electrodeChoice,arrayType,waveDetectionMethod);
        end
    end
end

%% get unique directions and circ correlations
overlapUqDir1 = cell(1,length(outputsTW1));
rhoValue = zeros(2,length(outputsTW1));
pValue = zeros(2,length(outputsTW1));
uqDirSG1 = cell(1,length(outputsTW1));
uqDirFG1 = cell(1,length(outputsTW1));

overlapUqDir2 = cell(1,length(outputsTW1));
uqDirSG2 = cell(1,length(outputsTW1));
uqDirFG2 = cell(1,length(outputsTW1));

for i = 1:length(outputsTW1)
    [dirSG,uniqueDirSG,sgBounds] = getWaveSegments(outputsTW{k,j},timeVals,wobbleLim,segmentationMethod,boundryLims,waveLengthLimit);
    [newBounds,allDirSg,allDirFg,allUniqueDirs] = getOverlappingWaves(dirSg,sgBounds,dirFg,fgBounds,overlap);
    overlapUqDir1{i} = cell2mat(allUniqueDirs);
    [rhoValue(1,i), pValue(1,i)] = circ_corrcc(overlapUqDir1{i}(1,:),overlapUqDir1{i}(2,:));
    allDirSg(isnan(allDirSg)) = [];
    uqDirSG1{i} = wrapToPi(unique(allDirSg));
    allDirFg(isnan(allDirFg)) = [];
    uqDirFG1{i} = wrapToPi(unique(allDirFg));
end
clear i allDirSg allDirFg allUniqueDirs 

for i = 1:numel(outputsTW2)
    [~,~,allDirSg,allDirFg,allUniqueDirs] = findOverlappingWaves(outputsTW2{i},timeVals,0.5,2);
    overlapUqDir2{i} = cell2mat(allUniqueDirs);
    [rhoValue(2,i), pValue(2,i)] = circ_corrcc(overlapUqDir2{i}(1,:),overlapUqDir2{i}(2,:));
    allDirSg(isnan(allDirSg)) = [];
    uqDirSG2{i} = wrapToPi(unique(allDirSg));
    allDirFg(isnan(allDirFg)) = [];
    uqDirFG2{i} = wrapToPi(unique(allDirFg));
end
clear i allDirSg allDirFg allUniqueDirs 

%% plot the data - polar plot with only the mean directions
tiledlayout(2,3,'TileSpacing','Compact');
colors = parula(numel(orientations));
colorVals = cat(1,[52 148 186]./255,[236 112 22]./255);

legends = {['0' char(176)], ['22.5' char(176)], ['45' char(176)], ['67.5' char(176)], ['90' char(176)], ['112.5' char(176)], ['135' char(176)], ['157.5' char(176)]};
meanReq = 2;
data1 = uqDirSG1;
data2 = uqDirFG1;
data3 = uqDirSG2;
data4 = uqDirFG2;

% plot directions
nexttile(1)
getViolinPlotPlain(data1,colors,legends,meanReq)
title('Slow Gamma Wave Directions')
ylabel('Directions')
% ylim([-3.14 3.14])
yticks([-3.14,0,3,14])
yticklabels({'-180','0','180'})

% violin plots monkey1 slow gamma - ori
nexttile(2)
getViolinPlotPlain(data2,colors,legends,meanReq)
title('Fast Gamma Wave Directions')
ylabel('Directions')
% ylim([-3.14 3.14])
yticks([-3.14,0,3,14])
yticklabels({'-180','0','180'})

% violin plots monkey2 slow gamma - ori
nexttile(4)
getViolinPlotPlain(data3,colors,legends,meanReq)
title('Slow Gamma Wave Directions')
ylabel('Directions')
% ylim([-3.14 3.14])
yticks([-3.14,0,3,14])
yticklabels({'-180','0','180'})
xlabel('Orientations')

% violin plots monkey2 slow gamma - ori
nexttile(5)
getViolinPlotPlain(data4,colors,legends,meanReq)
title('Fast Gamma Wave Directions')
ylabel('Directions')
% ylim([-3.14 3.14])
yticks([-3.14,0,3,14])
yticklabels({'-180','0','180'})
xlabel('Orientations')


%plot the corrcoef
nexttile(3)
sigPs = find(pValue(1,:)<0.05);
plot(1:numel(orientations),rhoValue(1,:),'-o','LineWidth',1.5,'MarkerIndices',1:numel(orientations),'Color',colorVals(1,:))
if ~isempty(sigPs)
    hold on
    sigPlot = ones(1,length(sigPs));
    sigPlot(sigPs) = 1;
%     sigPlot(sigPs1) = rhoValue(1,sigPs1)+rhoValue(1,sigPs1).*0.5;
    plot(1:numel(sigPlot),sigPlot*0.3,'*','Color','k','LineWidth',1.5)
end
ylim([-1 1])
xticks([1 2 3 4 5 6 7 8])
xticklabels(legends)
xlabel('Orientations')
ylabel('Circular CC')
title('Circ Corr across orientations:M1')

nexttile(6)
sigPs = find(pValue(2,:)<0.05);
plot(1:numel(orientations),rhoValue(2,:),'-o','LineWidth',1.5,'MarkerIndices',1:numel(orientations),'Color',colorVals(1,:))
if ~isempty(sigPs)
    hold on
    sigPlot = ones(1,length(sigPs));
    sigPlot(sigPs) = 1;
%     sigPlot(sigPs1) = rhoValue(1,sigPs1)+rhoValue(1,sigPs1).*0.5;
    plot(1:numel(sigPlot),sigPlot*0.3,'*','Color','k','LineWidth',1.5)
end
ylim([-1 1])
xticks([1 2 3 4 5 6 7 8])
xticklabels(legends)
xlabel('Orientations')
ylabel('Circular CC')
title('Circ Corr across orientations:M2')
