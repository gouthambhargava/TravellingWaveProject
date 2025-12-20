%generate fig 4 - final
% set up some parameters
outputsTW1 = {load('alpaH_12_0.5T_selected_met1.mat'),load('alpaH_22_0.5T_selected_met1.mat'),load('alpaH_32_0.5T_selected_met1.mat'),load('alpaH_42_0.5T_selected_met1.mat'),load('alpaH_52_0.5T_selected_met1.mat'),load('alpaH_62_0.5T_selected_met1.mat'),load('alpaH_72_0.5T_selected_met1.mat'),load('alpaH_82_0.5T_selected_met1.mat')};
outputsTW2 = {load('kesariH_12_0.5T_selected_met1.mat'),load('kesariH_22_0.5T_selected_met1.mat'),load('kesariH_32_0.5T_selected_met1.mat'),load('kesariH_42_0.5T_selected_met1.mat'),load('kesariH_52_0.5T_selected_met1.mat'),load('kesariH_62_0.5T_selected_met1.mat'),load('kesariH_72_0.5T_selected_met1.mat'),load('kesariH_82_0.5T_selected_met1.mat')};

thresh = 0.5;
electrodeFraction = 0.5;
electrodeChoice = 'selected';
arrayType = 'Microelectrode';
waveDetectionMethod = 1;
orientations = [0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5];
overlap = 0.5;

freqRangeList{1} = [20 35]; freqRangeList{2} = [40 60];
numFrequencyRanges = numel(freqRangeList);
sPos = 2; % spatial frequency: 0.5 (1), 1(2), 2 (3), 4 (4), 8 (5), all SFs (6). Note that the same code can be used for the size project also later where stimulus size is changed instead of spatial frequency
oriPos = 1:8; % orientation: 0 (1), 22.5 (2), 45 (3), 67.5 (4), 90 (5), 112.5 (6), 135 (7), 157.5 (8), all orientations (9)
stimPeriod = [0.25 0.75];

%% load data and get wave parameters
%for monkey 1
wobble1 = 5; % in deg
wobble2 = 10;
lengthLimit = 10;% in ms
segmentationMethod = 2;

% outputsTW1 = cell(1,length(oriPos));
rhoValM1 = nan(length(oriPos),2);
pValM1 = nan(length(oriPos),2);
rhoValM2 = nan(length(oriPos),2);
pValM2 = nan(length(oriPos),2);

for i = 1:length(oriPos) % for all ori
    [~,rhoValM1(i,1),pValM1(i,1)] = simpleOvWrapper(outputsTW1{i},timeVals,wobble1,segmentationMethod,lengthLimit,0.5);
    [~,rhoValM1(i,2),pValM1(i,2)] = simpleOvWrapper(outputsTW1{i},timeVals,wobble2,segmentationMethod,lengthLimit,0.5);
    [~,rhoValM2(i,1),pValM2(i,1)] = simpleOvWrapper(outputsTW2{i},timeVals,wobble1,segmentationMethod,lengthLimit,0.5);
    [~,rhoValM2(i,2),pValM2(i,2)] = simpleOvWrapper(outputsTW2{i},timeVals,wobble2,segmentationMethod,lengthLimit,0.5);
end
%% get unique directions and circ correlations
wobbleLim = 0;
segmentationMethod = 3;
waveLengthLimit = 10;

rhoValue = zeros(2,length(outputsTW1));
pValue = zeros(2,length(outputsTW1));
sgDirsM1 = cell(1,length(outputsTW1));
fgDirsM1 = cell(1,length(outputsTW1));

for i = 1:length(outputsTW1)
    outputs = outputsTW1{i};
    numTrials = length(outputs);
    waveVector = nan(numTrials,length(timeVals),numFrequencyRanges);
    waveBounds = cell(numFrequencyRanges,numTrials);
    for j = 1:numFrequencyRanges
        for k = 1:numTrials
            [waveVector(k,:,j),~,waveBounds{j,k}] = getWaveSegments(outputs{j,k},timeVals,wobbleLim,segmentationMethod,stimPeriod,waveLengthLimit);
        end
    end
    allDirSg = nan(numTrials,length(timeVals));
    allDirFg = nan(numTrials,length(timeVals));
    allUniqueDirs = cell(1,numTrials);
    for k = 1:numTrials
        [~,allDirSg(k,:),allDirFg(k,:),allUniqueDirs{k}] = getOverlappingWaves(waveVector(k,:,1),waveBounds{1,k},waveVector(k,:,2),waveBounds{2,k},overlap);
    end
    ccDirs = cell2mat(allUniqueDirs);
    ccDirs(:,isnan(ccDirs(1,:))) = [];
    [rhoValue(1,i), pValue(1,i)] = circ_corrcc(ccDirs(1,:),ccDirs(2,:));
    sgDirsM1{i} = allDirSg(~isnan(allDirSg));
    fgDirsM1{i} = allDirFg(~isnan(allDirFg));
end
clear i allDirSg allDirFg allUniqueDirs cc1 cc2 waveVector waveBounds i j k 


sgDirsM2 = cell(1,length(outputsTW1));
fgDirsM2 = cell(1,length(outputsTW1));


for i = 1:length(outputsTW2)
    outputs = outputsTW2{i};
    numTrials = length(outputs);
    waveVector = nan(numTrials,length(timeVals),numFrequencyRanges);
    waveBounds = cell(numFrequencyRanges,numTrials);
    for j = 1:numFrequencyRanges
        for k = 1:numTrials
            [waveVector(k,:,j),~,waveBounds{j,k}] = getWaveSegments(outputs{j,k},timeVals,wobbleLim,segmentationMethod,stimPeriod,waveLengthLimit);
        end
    end
    allDirSg = nan(numTrials,length(timeVals));
    allDirFg = nan(numTrials,length(timeVals));
    allUniqueDirs = cell(1,numTrials);
    for k = 1:numTrials
        [~,allDirSg(k,:),allDirFg(k,:),allUniqueDirs{k}] = getOverlappingWaves(waveVector(k,:,1),waveBounds{1,k},waveVector(k,:,2),waveBounds{2,k},overlap);
    end
    ccDirs = cell2mat(allUniqueDirs);
    ccDirs(:,isnan(ccDirs(1,:))) = [];
    [rhoValue(2,i), pValue(2,i)] = circ_corrcc(ccDirs(1,:),ccDirs(2,:));
    sgDirsM2{i} = allDirSg(~isnan(allDirSg));
    fgDirsM2{i} = allDirFg(~isnan(allDirFg));
end
clear i allDirSg allDirFg allUniqueDirs cc1 cc2 waveVector waveBounds i j k 


%% plot the data - polar plot with only the mean directions
tiledlayout(2,3,'TileSpacing','Compact');
colors = parula(numel(orientations));
colorVals = cat(1,[52 148 186]./255,[236 112 22]./255);

legends = {['0' char(176)], ['22.5' char(176)], ['45' char(176)], ['67.5' char(176)], ['90' char(176)], ['112.5' char(176)], ['135' char(176)], ['157.5' char(176)]};


% plot mean angles for slow gamma M1
nexttile(1)
getMultiMeanPlot(sgDirsM1,colors)
title('Slow Gamma Wave Directions')


% plot mean angles for fast gamma M1
nexttile(2)
getMultiMeanPlot(fgDirsM1,colors)
title('Fast Gamma Wave Directions')


% plot mean angles for slow gamma M2
nexttile(4)
getMultiMeanPlot(fgDirsM2,colors)
title('Slow Gamma Wave Directions')


% plot mean angles for fast gamma M2
nexttile(5)
getMultiMeanPlot(sgDirsM2,colors)
title('Fast Gamma Wave Directions')

pvalSig = 0.05/8;
% plot of circ correlations M1
nexttile(3)
sigPs = find(pValue(1,:)<pvalSig);
plot(1:numel(orientations),rhoValue(1,:),'-o','LineWidth',1.5,'MarkerIndices',1:numel(orientations),'Color',colorVals(1,:))
if ~isempty(sigPs)
    hold on
    sigPlot = ones(1,length(sigPs));
    sigPlot(sigPs) = 1;
    plot(1:numel(sigPlot),sigPlot*0.3,'*','Color',colorVals(1,:),'LineWidth',1.5)
end
hold on
sigPs = find(pValM1(:,1)<pvalSig);
plot(1:numel(orientations),rhoValM1(:,1),'-o','LineWidth',1.5,'MarkerIndices',1:numel(orientations),'Color','red')
if ~isempty(sigPs)
    hold on
    sigPlot = ones(1,length(sigPs));
    sigPlot(sigPs) = 1;
    plot(1:numel(sigPlot),sigPlot*0.3,'*','Color','red','LineWidth',1.5)
end

sigPs = find(pValM1(:,2)<pvalSig);
plot(1:numel(orientations),rhoValM1(:,2),'-o','LineWidth',1.5,'MarkerIndices',1:numel(orientations),'Color','green')
if ~isempty(sigPs)
    hold on
    sigPlot = ones(1,length(sigPs));
    sigPlot(sigPs) = 1;
    plot(1:numel(sigPlot),sigPlot*0.3,'*','Color','green','LineWidth',1.5)
end


ylim([-1 1])
xticks([1 2 3 4 5 6 7 8])
xticklabels(legends)
xlabel('Orientations')
ylabel('Circular CC')
title('Circ Corr across orientations:M1')

% plot of circ correlations M2
nexttile(6)
sigPs = find(pValue(2,:)<pvalSig);
plot(1:numel(orientations),rhoValue(2,:),'-o','LineWidth',1.5,'MarkerIndices',1:numel(orientations),'Color',colorVals(1,:))
if ~isempty(sigPs)
    hold on
    sigPlot = ones(1,length(sigPs));
    sigPlot(sigPs) = 1;
    plot(1:numel(sigPlot),sigPlot*0.3,'*','Color',colorVals(1,:),'LineWidth',1.5)
end
hold on
sigPs = find(pValM2(:,1)<pvalSig);
plot(1:numel(orientations),rhoValM2(:,1),'-o','LineWidth',1.5,'MarkerIndices',1:numel(orientations),'Color','red')
if ~isempty(sigPs)
    hold on
    sigPlot = ones(1,length(sigPs));
    sigPlot(sigPs) = 1;
    plot(1:numel(sigPlot),sigPlot*0.3,'*','Color','red','LineWidth',1.5)
end

sigPs = find(pValM2(:,2)<pvalSig);
plot(1:numel(orientations),rhoValM2(:,2),'-o','LineWidth',1.5,'MarkerIndices',1:numel(orientations),'Color','green')
if ~isempty(sigPs)
    hold on
    sigPlot = ones(1,length(sigPs));
    sigPlot(sigPs) = 1;
    plot(1:numel(sigPlot),sigPlot*0.3,'*','Color','green','LineWidth',1.5)
end

ylim([-1 1])
xticks([1 2 3 4 5 6 7 8])
xticklabels(legends)
xlabel('Orientations')
ylabel('Circular CC')
title('Circ Corr across orientations:M2')

annotation('textarrow',[0.36 0.38],[0.9 0.9],'String',legends{1},'Color',colors(1,:),'LineWidth',2)
annotation('textarrow',[0.36 0.38],[0.87 0.87],'String',legends{2},'Color',colors(2,:),'LineWidth',2)
annotation('textarrow',[0.36 0.38],[0.84 0.84],'String',legends{3},'Color',colors(3,:),'LineWidth',2)
annotation('textarrow',[0.36 0.38],[0.81 0.81],'String',legends{4},'Color',colors(4,:),'LineWidth',2)
annotation('textarrow',[0.36 0.38],[0.78 0.78],'String',legends{5},'Color',colors(5,:),'LineWidth',2)
annotation('textarrow',[0.36 0.38],[0.75 0.75],'String',legends{6},'Color',colors(6,:),'LineWidth',2)
annotation('textarrow',[0.36 0.38],[0.72 0.72],'String',legends{7},'Color',colors(7,:),'LineWidth',2)
annotation('textarrow',[0.36 0.38],[0.69 0.69],'String',legends{8},'Color',colors(8,:),'LineWidth',2)

annotation('textbox',[0.09, 0.95, 0.02, 0.02], 'string', 'A','FontSize',25,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.37, 0.95, 0.02, 0.02], 'string', 'B','FontSize',25,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.61, 0.95, 0.02, 0.02], 'string', 'C','FontSize',25,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.09, 0.5, 0.02, 0.02], 'string', 'D','FontSize',25,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.37, 0.5, 0.02, 0.02], 'string', 'E','FontSize',25,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.61, 0.5, 0.02, 0.02], 'string', 'F','FontSize',25,'FontWeight','bold','EdgeColor','none')

annotation('textbox',[0.93, 0.82, 0.02, 0.02], 'string', ['0',char(176)],'FontSize',15,'FontWeight','bold','EdgeColor','none','Color',colorVals(1,:))
annotation('textbox',[0.93, 0.77, 0.02, 0.02], 'string', ['5',char(176)],'FontSize',15,'FontWeight','bold','EdgeColor','none','Color','red')
annotation('textbox',[0.93, 0.72, 0.02, 0.02], 'string', ['10',char(176)],'FontSize',15,'FontWeight','bold','EdgeColor','none','Color','green')

annotation('textbox',[0.93, 0.34, 0.02, 0.02], 'string', ['0',char(176)],'FontSize',15,'FontWeight','bold','EdgeColor','none','Color',colorVals(1,:))
annotation('textbox',[0.93, 0.29, 0.02, 0.02], 'string', ['5',char(176)],'FontSize',15,'FontWeight','bold','EdgeColor','none','Color','red')
annotation('textbox',[0.93, 0.24, 0.02, 0.02], 'string', ['10',char(176)],'FontSize',15,'FontWeight','bold','EdgeColor','none','Color','green')
