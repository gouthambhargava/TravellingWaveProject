%% get video TW
%% for alpa
sPos = 2;
oriPos = 4;
req = 1;
dataPath = 'F:\monkeyData\data';
gridType = 'Microelectrode';
subjectName='alpaH'; expDate = '210817'; protocolName = 'GRF_002';
[mData,goodElectrodes,timeVals,~,~] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sPos,oriPos);
load('alpaH_42_0.5T_selected_met1.mat');
outputsTW = outputs;
%%
numTrials = size(mData,2);
freqRangeList{1} = [20 35]; freqRangeList{2} = [40 60];
numGoodElectrodes = length(goodElectrodes);
numFrequencyRanges = numel(freqRangeList);


% Do burst estimation
thresholdFactor = 3;
% stimulusDurationS = [0 0.8]; % Stimulus duration to be highlighted
baselinePeriodS = [-0.5 0];
stimulusPeriodS = [0.25 0.75];
analysisPeriodS = [-0.5 1];
filterOrder = 4;
for iFreq=1:numFrequencyRanges
    for iElec=1:numGoodElectrodes
        [~,~,~,burstTS(iElec,:,:,iFreq),filteredSignal(iElec,:,:,iFreq),~] = getHilbertBurst(squeeze(mData(iElec,:,:)),timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,freqRangeList{iFreq},filterOrder,1,analysisPeriodS);
    end
end

%%
trial = 22;
wave = 3;
wobbleLim = 0;
segOption = 2;
boundryLims = [0.25 0.75];
lengthLimit = 10;

for i = 1:numTrials
    for j = 1:numFrequencyRanges
        [waveVector(i,:,j),~,waveBounds{j,i}] = getWaveSegments(outputsTW{j,i},timeVals,wobbleLim,segOption,boundryLims, lengthLimit);
        directions(:,:,i,j) = outputsTW{j,i}.direction;
        phases(:,:,i,j) = angle(hilbert(squeeze(filteredSignal(:,i,:,j))'))';
    end
end

% find overlapping waves
allDirSG = nan(numTrials,length(timeVals));
allDirFG = nan(numTrials,length(timeVals));
newBounds = cell(1,numTrials);
uniqueDirs = cell(1,numTrials);
overlap = 0.5;
for i = 1:numTrials
    [newBounds{i},allDirSG(i,:),allDirFG(i,:),uniqueDirs{i},~,intPts(i,:)] = getOverlappingWaves(waveVector(i,:,1),waveBounds{1,i},waveVector(i,:,2),waveBounds{2,i},overlap);
end


pgd1 = outputsTW{1,trial}.pgd(1,:);
pgd2 = outputsTW{2,trial}.pgd(1,:);
pgd1(pgd1<0) = 0;
pgd1(isnan(pgd1)) = 0;
pgd2(pgd2<0) = 0;
pgd2(isnan(pgd2)) = 0;

duration1 = newBounds{1,trial}{1,1}(:,wave);
duration2 = newBounds{1,trial}{2,1}(:,wave);
durOverlap = intersect(duration1(1):duration1(2),duration2(1):duration2(2));


gridLayout = rot90(reshape(1:81,[9,9]),2); %set the grid layout
phaseData1 = nan(size(gridLayout,1),size(gridLayout,2),numel(durOverlap));
phaseData2 = nan(size(gridLayout,1),size(gridLayout,2),numel(durOverlap));
dirData1 = nan(size(gridLayout,1),size(gridLayout,2),numel(durOverlap));
dirData2 = nan(size(gridLayout,1),size(gridLayout,2),numel(durOverlap));
burstFrac = squeeze(nansum(squeeze(burstTS(:,trial,:,:)),1));
burstFrac(burstFrac<0.5*numGoodElectrodes) = nan;
burstFrac = burstFrac/numGoodElectrodes;

for i = 1:numel(goodElectrodes)
    [x,y] = find(gridLayout==goodElectrodes(i));
    phaseData1(x,y,:) = angle(hilbert(squeeze(filteredSignal(i,trial,durOverlap,1))));
    phaseData2(x,y,:) = angle(hilbert(squeeze(filteredSignal(i,trial,durOverlap,2))));
    dirData1(x,y,:) = allDirSG(trial,durOverlap);
    dirData2(x,y,:) = allDirFG(trial,durOverlap);
end


timeValues = timeVals(durOverlap(1):durOverlap(end));
u1 = cos(dirData1);
v1 = sin(dirData1);
u2 = cos(dirData2);
v2 = sin(dirData2);
%%
f = figure;
set(gcf, 'Position', get(0, 'Screensize'));
% f.Position = [100 100 540 400];
colorVals = cat(1,[52 148 186]./255,[236 112  22]./255);
[x, y] = meshgrid(1:9,1:9);
pgd1(isnan(pgd1)) = 0;
pgd2(isnan(pgd2)) = 0;

a = waveVector(trial,:,1);
b = waveVector(trial,:,2);
c = intPts(trial,:);
c(~isnan(c)) = 1;
a(~isnan(a)) = 1;
b(~isnan(b)) = 1;

pgdGrid = getPlotHandles(1,1,[0.1 0.7 0.8 0.25],0.01,0.01,0);
phaseGrid1 = getPlotHandles(1,1,[0.05 0.05 0.4 0.5],0.01,0.01,0);
phaseGrid2 = getPlotHandles(1,1,[0.55 0.05 0.4 0.5],0.01,0.01,0);

plot(pgdGrid,timeVals,pgd1,'LineWidth',1,'Color',colorVals(1,:))
hold(pgdGrid,"on")
plot(pgdGrid,timeVals,pgd2,'LineWidth',1,'Color',colorVals(2,:))
plot(pgdGrid,timeVals,a*-0.1,'|','Color',colorVals(1,:))
plot(pgdGrid,timeVals,b*-0.2,'|','Color',colorVals(2,:))
plot(pgdGrid,timeVals,c*-0.3,'|','Color','black')
title(pgdGrid,'PGD of overlapping TWs:M1')

frames = [];
writerObj = VideoWriter('SuppVideo1.avi');
writerObj.FrameRate = 10;
% open the video writer
open(writerObj);
for i = 1:numel(durOverlap)
    h = line([timeValues(i),timeValues(i)],[-0.4,1],'Parent',pgdGrid,'Color','black');
    legend(pgdGrid,'Slow Gamma','Fast Gamma','location','best')
    xlim(pgdGrid,[0 0.8])
    imagesc(cos(phaseData1(:,:,i)),'Parent',phaseGrid1)
    clim(phaseGrid1,[-1 1])
    hold(phaseGrid1,'on')
    quiver(x,y,u1(:,:,i),v1(:,:,i),'color','white','AutoScaleFactor',0.9,'parent',phaseGrid1)
    set(phaseGrid1, 'YDir', 'Normal')
    axis(phaseGrid1,'square')
    axis(phaseGrid1, 'off')
    title(phaseGrid1,['Slow Gamma:', num2str(timeValues(i)),'s'])

    imagesc(cos(phaseData2(:,:,i)),'Parent',phaseGrid2)
    clim(phaseGrid2,[-1 1])
    hold(phaseGrid2,'on')
    quiver(x,y,u2(:,:,i),v2(:,:,i),'color','white','AutoScaleFactor',0.9,'parent',phaseGrid2)
    set(phaseGrid2, 'YDir', 'Normal')
    axis(phaseGrid2,'square')
    axis(phaseGrid2, 'off')
    title(phaseGrid2,['Fast Gamma:', num2str(timeValues(i)),'s'])
    % drawnow
    pause(0.1)
    frames = getframe(gcf);
    writeVideo(writerObj, frames);

    % exportgraphics(gcf,'Supplementary Video 1.gif','Append',true);
    delete(h)
end
