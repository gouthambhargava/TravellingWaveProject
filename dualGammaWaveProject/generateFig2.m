%% generate fig 2
% load data
dataPath = 'G:\monkeyData\data';
gridType = 'Microelectrode';
subjectName='alpaH'; expDate = '210817'; protocolName = 'GRF_002'; 
sPos = 2; % spatial frequency: 0.5 (1), 1(2), 2 (3), 4 (4), 8 (5), all SFs (6). Note that the same code can be used for the size project also later where stimulus size is changed instead of spatial frequency
oriPos = 4; % orientation: 0 (1), 22.5 (2), 45 (3), 67.5 (4), 90 (5), 112.5 (6), 135 (7), 157.5 (8), all orientations (9)
[mData,goodElectrodes,timeVals,rfData,parameters] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sPos,oriPos);
freqRangeList{1} = [20 35]; freqRangeList{2} = [40 60];
% get burst data 
% burst calculation parameters
thresholdFactor = 3;
baselinePeriodS = [-0.5 0];
stimulusPeriodS = [0.25 0.75];
analysisPeriodS = [-0.5 1];
filterOrder = 4;
numFrequencyRanges = numel(freqRangeList);

burstTS = zeros(numel(goodElectrodes),size(mData,2),size(mData,3),numFrequencyRanges);
filteredSignal = zeros(numel(goodElectrodes),size(mData,2),size(mData,3),numFrequencyRanges);
for iFreq=1:numFrequencyRanges
    for iElec=1:numel(goodElectrodes)
        [~,~,~,burstTS(iElec,:,:,iFreq),filteredSignal(iElec,:,:,iFreq)] = getHilbertBurst(squeeze(mData(iElec,:,:)),timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,freqRangeList{iFreq},filterOrder,1,analysisPeriodS);
     end
end
bursts = burstTS;
bursts(isnan(bursts)) = 0;
burstFrac = squeeze(sum(bursts))*100/numel(goodElectrodes);
burstFrac(burstFrac<thresh) = 0;
burstFrac(burstFrac==0) = nan;
burstFrac(~isnan(burstFrac)) = 1;

% get outputs for M1
% get list of electrode locations
numGoodElectrodes = numel(goodElectrodes);
electrodeArray = rot90(reshape(1:81,[9,9]),2); %set the grid layout
locList = nan(numGoodElectrodes,2);
for iElec = 1:numGoodElectrodes
    [locList(iElec,1),locList(iElec,2)] = find(electrodeArray==goodElectrodes(iElec));
end

% calculate TW parameters
tic
outputsTW = cell(numFrequencyRanges,numTrials);
for i = 1:numTrials
    for j = 1:numFrequencyRanges
        burstMat = squeeze(burstTS(:,i,:,j));
        phiMat = angle(hilbert(squeeze(filteredSignal(:,i,:,j))'))';
        outputsTW{j,i} = getTWCircParams(phiMat,burstMat,timeVals,goodElectrodes,locList,electrodeFraction,'selected');
    end
end
toc
% define some parameters for wave detection
lengthLimit = 25; %ms
boundryLims = [0.25 0.75];
wobbleLim = 5; %degree
segOption = 3;
numFrequencyRanges = size(outputsTW,1);

% segment waves and get overlapping waves for all trials
%initialize outputs
waveVector = nan(numTrials,length(timeVals),numFrequencyRanges);
uniqueDirs = cell(numFrequencyRanges,numTrials);
waveBounds = cell(numFrequencyRanges,numTrials);
for i = 1:numTrials
    for j = 1:numFrequencyRanges
    [waveVector(i,:,j),uniqueDirs{j,i},waveBounds{j,i}] = getWaveSegments(outputsM1{j,i},timeVals,wobbleLim,segOption,boundryLims, lengthLimit);
    end
end

% find overlapping waves 
allUniqueDirs = [];
dirSG = nan(numTrials,length(timeVals));
dirFG = nan(numTrials,length(timeVals));
overlap = 0.5;
for i = 1:numTrials
    [~,dirSG(i,:),dirFG(i,:),uniqueDirs] = getOverlappingWaves(waveVector(i,:,1),waveBounds{1,i},waveVector(i,:,2),waveBounds{2,i},overlap);
    allUniqueDirs = cat(2,allUniqueDirs,uniqueDirs);
end

%% get overlapping TW's
wobble = 5;
trial = 15;
wave = 1;

[newBounds,sigCells,allDirSg,allDirFg,uqDir] = findOverlappingWaves(outputsTW,timeVals,0.5,wobble);
[sgDir,~,boundsAllSG] = getWaveSegments(outputsTW,timeVals,wobble);
[fgDir,~,boundsAllFG] = getWaveSegments(outputsTW,timeVals,wobble);
overlappingDirsSg = nan(size(sgDir));
overlappingDirsFg = nan(size(fgDir));
overlappingDirsSg(sigCells,:) = allDirSg;
overlappingDirsFg(sigCells,:) = allDirFg;

pgd(:,1) = outputsTW{1,sigCells(trial)}.pgd;
pgd(:,2) = outputsTW{2,sigCells(trial)}.pgd;
pgd(pgd<0) = 0;
pgd(isnan(pgd)) = 0;

directionO(:,1) = allDirSg(trial,:);
directionO(:,2) = allDirFg(trial,:);
direction(:,1) = sgDir(sigCells(trial),:);
direction(:,2) = fgDir(sigCells(trial),:);

direction(~isnan(directionO)) = nan;

durIndices = [newBounds{1,trial}(:,wave),newBounds{2,trial}(:,wave)];
durIndices = intersect(durIndices(1,1):durIndices(2,1),durIndices(1,2):durIndices(2,2));

selectedDir{1} = directionO(durIndices,1);
selectedDir{2} = directionO(durIndices,2);

% generate phase data
% get phase data for plot 3
gridLayout = rot90(reshape(1:81,[9,9]),2); %set the grid layout
phaseData1 = nan(size(gridLayout,1),size(gridLayout,2),numel(durIndices));
phaseData2 = nan(size(gridLayout,1),size(gridLayout,2),numel(durIndices));

for i = 1:numel(goodElectrodes)
    [x,y] = find(gridLayout==goodElectrodes(i));
    phaseData1(x,y,:) = angle(hilbert(squeeze(filteredSignal(i,sigCells(trial),durIndices,1))));
    phaseData2(x,y,:) = angle(hilbert(squeeze(filteredSignal(i,sigCells(trial),durIndices,2))));
end

[x, y] = meshgrid(1:9,1:9);
u1 = cos(selectedDir{1});
v1 = sin(selectedDir{1});
u2 = cos(selectedDir{2});
v2 = sin(selectedDir{2});

phaseInd1 = floor(linspace(48,85,6));
phaseInd2 = floor(linspace(48,85,6));
timeValsInt = timeVals(durIndices);
timeValsInt = [timeValsInt(phaseInd1(1)),timeValsInt(phaseInd1(end))];

%% get background gradient with imagesc for TW progression - for figure 2 - slow gamma
figure(1)
xVals = [0 100 100 53 52.5 0];
yVals = [0 0 100 120 120 100];
vertices = cat(1,xVals,yVals)';
faces = 1:numel(xVals);
colorValsNo = cat(1,[52 148 186]./255,[236 112 22]./255);
colorVals = cat(1,[0.1412,0,0.9412],[0.8 0 0]);
factors = linspace(0,1,numel(xVals));
for i = 1:numel(factors)
    colors(i,:) = colorVals(1,:)+(1-colorVals(1,:))*factors(i);
end
% colors = colorVals(1,:).*factors';

colors = flipud(colors);
patch('Vertices',vertices,'Faces',faces,'FaceVertexCData',colors,'FaceColor','interp','EdgeColor','none','FaceAlpha',.7)

axis off
ylim([0 150])
% plot pgd and burstFrac
axes('Position',[0.13 0.763 0.775 0.15])
burstFracRed = squeeze(burstFrac(sigCells(trial),:,:));
burstFracRed(burstFracRed==0) = nan;
burstFracRed(~isnan(burstFracRed)) = 1;
% x1 = [timeVals(durIndices(1,1)) timeVals(durIndices(2,1)) timeVals(durIndices(2,1)) timeVals(durIndices(1,1))];
y1 = [-1 -1 1 1];
x1 = [timeValsInt(1),timeValsInt(end),timeValsInt(end),timeValsInt(1)];
patch(x1, y1, colorVals(1,:),'FaceAlpha',0.3,'EdgeColor',colorVals(1,:),'LineStyle','none')
hold on
plot(timeVals,pgd(:,1),'color','k','linewidth',1);
hold on
plot(timeVals,burstFracRed(:,1)*-0.2,'Marker','|','LineWidth',1.2,'Color',[0.4660 0.6740 0.1880])
temp = direction(:,1);
temp(~isnan(temp)) = 1;
plot(timeVals,temp*-0.05,'|','LineWidth',1.2,'Color',colorValsNo(1,:))
temp = directionO(:,1);
temp(~isnan(temp)) = 1;
plot(timeVals,temp*-0.05,'|','LineWidth',1.2,'Color',colorVals(1,:))


ylabel('PGD')
xlabel('Time(s)')
xlim([0.25 0.75])
ylim([-0.2 0.7])
title('Slow Gamma')

% plot phase progression
posD = [0.430 0.430 0.430 0.14 0.14 0.14];
posL = repmat([0.18 0.42 0.65],1,2);
sizePos = 0.2;
dataLabels = timeVals(durIndices(1):durIndices(end));
for i = 1:length(posD)
    axes('Position',[posL(i) posD(i) sizePos sizePos])
    imagesc(cos(phaseData1(:,:,phaseInd1(1,i))))
    hold on
    caxis([-1 1])
    quiver(x,y,ones(size(x)).*v1(phaseInd1(i)),ones(size(y)).*u1(phaseInd1(i)),'color','white','AutoScaleFactor',0.9)
    title(['Time:',num2str(round(dataLabels(phaseInd1(1,i)),3)),'s'])
    axis off
end
annotation('textbox',[0.08,0.95, 0, 0], 'string', 'B','FontSize',30,'FontWeight','bold')
%% get background gradient with imagesc for TW progression - for figure 2 - fast gamma
figure(2)
xVals = [0 100 100 53 52.5 0];
yVals = [0 0 100 120 120 100];
vertices = cat(1,xVals,yVals)';
faces = 1:numel(xVals);
% colorVals = cat(1,[52 148 186]./255,[236 112 22]./255);

factors = linspace(0,1,numel(xVals));
for i = 1:numel(factors)
    colors(i,:) = colorVals(2,:)+(1-colorVals(2,:))*factors(i);
end
% colors = colorVals(1,:).*factors';
colors = flipud(colors);
patch('Vertices',vertices,'Faces',faces,'FaceVertexCData',colors,'FaceColor','interp','EdgeColor','none','FaceAlpha',.7)
axis off
ylim([0 150])
% plot pgd and burstFrac
axes('Position',[0.13 0.763 0.775 0.15])
% x1 = [timeVals(durIndices(1,2)) timeVals(durIndices(2,2)) timeVals(durIndices(2,2)) timeVals(durIndices(1,2))];
x1 = [timeValsInt(1),timeValsInt(end),timeValsInt(end),timeValsInt(1)];

y1 = [-1 -1 1 1];
patch(x1, y1, colorVals(2,:),'FaceAlpha',0.3,'EdgeColor',colorVals(2,:),'LineStyle','none')
hold on
plot(timeVals,pgd(:,2),'color','k','linewidth',1);
hold on
plot(timeVals,burstFracRed(:,2)*-0.2,'|','LineWidth',1.2,'Color','green')
temp = direction(:,2);
temp(~isnan(temp)) = 1;
plot(timeVals,temp*-0.05,'|','LineWidth',1.2,'Color',colorValsNo(2,:))
temp = directionO(:,2);
temp(~isnan(temp)) = 1;
plot(timeVals,temp*-0.05,'|','LineWidth',1.2,'Color',colorVals(2,:))

ylabel('PGD')
xlabel('Time(s)')
xlim([0.25 0.75])
ylim([-0.2 0.9])
title('Fast Gamma')

% plot phase progression
posD = [0.430 0.430 0.430 0.14 0.14 0.14];
posL = repmat([0.18 0.42 0.65],1,2);
sizePos = 0.2;
dataLabels = timeVals(durIndices(1):durIndices(end));
for i = 1:length(posD)
    axes('Position',[posL(i) posD(i) sizePos sizePos])
    imagesc(cos(phaseData2(:,:,phaseInd2(1,i))))
    hold on
    clim([-1 1])
    quiver(x,y,ones(size(x)).*v2(phaseInd2(1,i)),ones(size(y)).*u2(phaseInd2(1,i)),'color','white','AutoScaleFactor',0.7)
    title(['Time:',num2str(round(dataLabels(phaseInd2(1,i)),3)),'s'])
    axis off
end

annotation('textbox',[0.08,0.95, 0, 0], 'string', 'C','FontSize',30,'FontWeight','bold')

%% generate fig 2 top panel
colorValsNo = cat(1,[52 148 186]./255,[236 112 22]./255);
colorVals = cat(1,[0.1412,0,0.9412],[0.8 0 0]);

sgDir = sgDir(sigCells,:);
fgDir = fgDir(sigCells,:);
boundsAllSG = boundsAllSG(sigCells);
boundsAllFG = boundsAllFG(sigCells);

meanBoundsSG = cellfun(@mean,boundsAllSG,'UniformOutput',false);
meanBoundsSG = cellfun(@floor,meanBoundsSG,'UniformOutput',false);
meanBoundsFG = cellfun(@mean,boundsAllFG,'UniformOutput',false);
meanBoundsFG = cellfun(@floor,meanBoundsFG,'UniformOutput',false);

meanBoundsSGO = cellfun(@mean,newBounds(1,:),'UniformOutput',false);
meanBoundsSGO = cellfun(@floor,meanBoundsSGO,'UniformOutput',false);
meanBoundsFGO = cellfun(@mean,newBounds(2,:),'UniformOutput',false);
meanBoundsFGO = cellfun(@floor,meanBoundsFGO,'UniformOutput',false);

% sgDir(~isnan(overlappingDirsSg)) = nan;
% fgDir(~isnan(overlappingDirsFg)) = nan;
sgDir(~isnan(sgDir)) = 1;
fgDir(~isnan(fgDir)) = 1;
allDirSg(~isnan(allDirSg)) = 1;
allDirFg(~isnan(allDirFg)) = 1;


% numOverlap = numel(newBounds(1,:));
values1 = 1:1:size(sgDir,1);
values2 = values1-0.3;
values3 = values1+0.4;

%% apply patch to the selected trial
axTp = getPlotHandles(1,1,[0.1 0.1 0.3 0.86],0.1,0.1,0);

plot(timeVals,sgDir.*values1','Color',colorValsNo(1,:),'LineWidth',1.2)
hold on
plot(timeVals,allDirSg.*values1','Color',colorVals(1,:),'LineWidth',1.2)
for i = 1:size(sgDir,1)
    int = meanBoundsSG{i};
    plot(timeVals(int),ones(1,length(int)).*values1(i),'*','Color',colorValsNo(1,:),'LineWidth',1.2)
end
for i = 1:size(sgDir,1)
    int = meanBoundsSGO{i};
    plot(timeVals(int),ones(1,length(int)).*values1(i),'*','Color',colorVals(1,:),'LineWidth',1.2)
end


plot(timeVals,fgDir.*values2','Color',colorValsNo(2,:),'LineWidth',1.2)
plot(timeVals,allDirFg.*values2','Color',colorVals(2,:),'LineWidth',1.2)
for i = 1:size(fgDir,1)
    int = meanBoundsFG{i};
    plot(timeVals(int),ones(1,length(int)).*values2(i),'*','Color',colorValsNo(2,:),'LineWidth',1.2)
end
for i = 1:size(fgDir,1)
    int = meanBoundsFGO{i};
    plot(timeVals(int),ones(1,length(int)).*values2(i),'*','Color',colorVals(2,:),'LineWidth',1.2)
end

yline(values3,':','LineWidth',0.75)
% patch(xValues, yValues, [0.6350 0.0780 0.1840],'FaceAlpha',0.2)
ylim([0 20.5])
xlim([0.25 0.75])
yticks([1 5 10 15 20])
yticklabels({'1','5','10','15','20'});
xlabel('Time(s)')
ylabel('Trials')
title('Overlapping Travelling Waves across all trials')
xValues = [0.25 0.75 0.75 0.25];
yValues = [trial-0.5 trial-0.5 trial+0.5 trial+0.5];
patch(xValues, yValues, [0.6350 0.0780 0.1840],'FaceAlpha',0.2)
%%
annotation('rectangle',[0.61, 0.85, 0.03, 0.02],'FaceColor',colorValsNo(1,:),'EdgeColor','none')
annotation('rectangle',[0.61, 0.82, 0.03, 0.02],'FaceColor',colorValsNo(2,:),'EdgeColor','none')
annotation('rectangle',[0.61, 0.79, 0.03, 0.02],'FaceColor',colorVals(1,:),'EdgeColor','none')
annotation('rectangle',[0.61, 0.76, 0.03, 0.02],'FaceColor',colorVals(2,:),'EdgeColor','none')

annotation('textbox',[0.64, 0.787, 0.1, 0.1],'String','Slow Gamma','FontSize',20,'EdgeColor','none')
annotation('textbox',[0.64, 0.757, 0.1, 0.1],'String','Fast Gamma','FontSize',20,'EdgeColor','none')
annotation('textbox',[0.64, 0.727, 0.1, 0.1],'String','Slow Gamma Overlap','FontSize',20,'EdgeColor','none')
annotation('textbox',[0.64, 0.697, 0.1, 0.1],'String','Fast Gamma Overlap','FontSize',20,'EdgeColor','none')



annotation('textbox',[0.05,0.95, 0, 0], 'string', 'A','FontSize',30,'FontWeight','bold')