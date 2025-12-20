%% generate fig 2
% load data
% load('D:\IISC_work\TWGitScripts\TravellingWaveProject\dualGammaWaveProject\data\alpaHM1.mat')
% load('D:\IISC_work\gitScripts\FinalWobbleScripts\alpa24Data.mat')

% load data
dataPath = 'F:\monkeyData\data';
gridType = 'Microelectrode';
subjectName='alpaH'; expDate = '210817'; protocolName = 'GRF_002'; 
sPos = 2; % spatial frequency: 0.5 (1), 1(2), 2 (3), 4 (4), 8 (5), all SFs (6). Note that the same code can be used for the size project also later where stimulus size is changed instead of spatial frequency
oriPos = 4; % orientation: 0 (1), 22.5 (2), 45 (3), 67.5 (4), 90 (5), 112.5 (6), 135 (7), 157.5 (8), all orientations (9)
[mData,goodElectrodes,timeVals,rfData,parameters] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sPos,oriPos);
% mData = allData;
load('alpaH_42_0.5T_selected_met1.mat');
outputsTW = outputs;

thresh = 0.5;
freqRangeList{1} = [20 35]; freqRangeList{2} = [40 60];
% get burst data 
% burst calculation parameters
thresholdFactor = 3;
baselinePeriodS = [-0.5 0];
stimulusPeriodS = [0.25 0.75];
analysisPeriodS = [-0.5 1];
filterOrder = 4;
numFrequencyRanges = numel(freqRangeList);
numTrials = size(mData,2);

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

%% define some parameters for wave detection
trial = 22;
wave = 3;
lengthLimit = 10; %ms
boundryLims = [0.25 0.75];
wobbleLim = 0; %degree
segOption = 3;
numFrequencyRanges = size(outputsTW,1);

% segment waves and get overlapping waves for all trials
%initialize outputs
waveVector = nan(numTrials,length(timeVals),numFrequencyRanges);
% uniqueDirs = cell(numFrequencyRanges,numTrials);
waveBounds = cell(numFrequencyRanges,numTrials);
directions = nan(numGoodElectrodes,length(timeVals),numTrials,numFrequencyRanges);
phases = nan(numGoodElectrodes,length(timeVals),numTrials,numFrequencyRanges);

for i = 1:numTrials
    for j = 1:numFrequencyRanges
        [waveVector(i,:,j),~,waveBounds{j,i}] = getWaveSegments(outputsTW{j,i},timeVals,wobbleLim,segOption,boundryLims, lengthLimit);
        directions(:,:,i,j) = outputsTW{j,i}.direction;
        phases(:,:,i,j) = angle(hilbert(squeeze(filteredSignal(:,i,:,j))'))';
    end
end

% find overlapping waves 
uniqueDirs = cell(1,numTrials);
dirSG = nan(numTrials,length(timeVals));
dirFG = nan(numTrials,length(timeVals));
overlap = 0.5;
for i = 1:numTrials
    [waveBoundsOv{i},dirSG(i,:),dirFG(i,:),uniqueDirs{1,i}] = getOverlappingWaves(waveVector(i,:,1),waveBounds{1,i},waveVector(i,:,2),waveBounds{2,i},overlap);
    % allUniqueDirs = cat(2,allUniqueDirs,uniqueDirs);
end

% plot all trials with overlapping waves
ovWaves = dirFG;
ovWaves(isnan(dirSG)) = nan;
ovWaves(~isnan(ovWaves)) = 1;

        
trialsNum = 0.6:numTrials;
trialLabels = 1:numTrials;
%
figure(1)
subplot(1,2,1)
plotWavesAllTrials(waveVector,timeVals,[0.25 0.75])
hold on
plot(timeVals,ovWaves.*trialsNum','LineWidth',0.8,'Color','black')
hold on
patch([0.25 0.75 0.75 0.25], [trial-1 trial-1 trial trial], 'red','FaceAlpha',0.2,'EdgeColor','red','LineStyle','none')
ylim([0 35])
ylabel('Trials')
xlabel('Time (s)')
yticks(0.5:1:numTrials)
yticklabels(num2str(trialLabels'))
title('Overlapping waves across all trials')
colorVals = cat(1,[52 148 186]./255,[236 112  22]./255);

annotation('textbox',[0.05,0.95, 0, 0], 'string', 'A','FontSize',30,'FontWeight','bold')
annotation('textbox',[0.01, 0.587, 0.1, 0.1],'String','Slow Gamma','FontSize',12,'Color',colorVals(1,:),'EdgeColor','none')
annotation('textbox',[0.01, 0.550, 0.1, 0.1],'String','Fast Gamma','FontSize',12,'Color',colorVals(2,:),'EdgeColor','none')
annotation('textbox',[0.01, 0.510, 0.1, 0.1],'String','Overlapping Waves','FontSize',12,'Color','black','EdgeColor','none')

% get overlapping TW's


pgd(:,1) = outputsTW{1,trial}.pgd(1,:);
pgd(:,2) = outputsTW{2,trial}.pgd(1,:);
pgd(pgd<0) = 0;
pgd(isnan(pgd)) = 0;

directionO(:,1) = dirSG(trial,:);
directionO(:,2) = dirFG(trial,:);

directionOv = directionO(:,1);
directionOv(isnan(directionO(:,2))) = nan;

durIndices = [waveBoundsOv{1,trial}{1,1}(:,wave),waveBoundsOv{1,trial}{2,1}(:,wave)];
durIndices = intersect(durIndices(1,1):durIndices(2,1),durIndices(1,2):durIndices(2,2));

selectedDir{1} = directionO(durIndices,1);
selectedDir{2} = directionO(durIndices,2);

% generate phase data
% get phase data for plot 3
gridLayout = rot90(reshape(1:81,[9,9]),2); %set the grid layout
phaseData1 = nan(size(gridLayout,1),size(gridLayout,2),numel(durIndices));
phaseData2 = nan(size(gridLayout,1),size(gridLayout,2),numel(durIndices));
dirData1 = nan(size(gridLayout,1),size(gridLayout,2),numel(durIndices));
dirData2 = nan(size(gridLayout,1),size(gridLayout,2),numel(durIndices));


for i = 1:numel(goodElectrodes)
    [x,y] = find(gridLayout==goodElectrodes(i));
    phaseData1(x,y,:) = angle(hilbert(squeeze(filteredSignal(i,trial,durIndices,1))));
    phaseData2(x,y,:) = angle(hilbert(squeeze(filteredSignal(i,trial,durIndices,2))));
    dirData1(x,y,:) = selectedDir{1};
    dirData2(x,y,:) = selectedDir{2};
end

[x, y] = meshgrid(1:9,1:9);

phaseInd1 = floor(linspace(1,length(durIndices),6));
phaseInd2 = floor(linspace(1,length(durIndices),6));
timeValsInt = timeVals(durIndices);
timeValsInt = [timeValsInt(1),timeValsInt(end)];

% dirData1 = dirData1+deg2rad(45);
% dirData2 = dirData2+deg2rad(100);
%% get background gradient with imagesc for TW progression - for figure 2 - slow gamma
figure(2)
xVals = [0 100 100 87 83 0];
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
burstFracRed = squeeze(burstFrac(trial,:,:));
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
temp = directionO(:,1);
temp(~isnan(temp)) = 1;
plot(timeVals,temp*-0.05,'|','LineWidth',1.2,'Color',colorValsNo(1,:))
temp = directionOv;
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
    clim([-1 1])
    quiver(x,y,cos(dirData1(:,:,i)),sin(dirData1(:,:,i)),'color','white','AutoScaleFactor',0.9)
    title(['Time:',num2str(round(dataLabels(phaseInd1(1,i)),3)),'s'])
    axis off
end
annotation('textbox',[0.08,0.95, 0, 0], 'string', 'B','FontSize',30,'FontWeight','bold')
% get background gradient with imagesc for TW progression - for figure 2 - fast gamma
figure(3)
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
temp = directionO(:,1);
temp(~isnan(temp)) = 1;
plot(timeVals,temp*-0.05,'|','LineWidth',1.2,'Color',colorValsNo(2,:))
temp = directionOv;
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
    quiver(x,y,cos(dirData2(:,:,i)),sin(dirData2(:,:,i)),'color','white','AutoScaleFactor',0.9)
    title(['Time:',num2str(round(dataLabels(phaseInd2(1,i)),3)),'s'])
    axis off
end

annotation('textbox',[0.08,0.95, 0, 0], 'string', 'C','FontSize',30,'FontWeight','bold')

%% generate fig 2 top panel

annotation('rectangle',[0.61, 0.85, 0.03, 0.02],'FaceColor',colorValsNo(1,:),'EdgeColor','none')
annotation('rectangle',[0.61, 0.82, 0.03, 0.02],'FaceColor',colorValsNo(2,:),'EdgeColor','none')
annotation('rectangle',[0.61, 0.79, 0.03, 0.02],'FaceColor',colorVals(1,:),'EdgeColor','none')
annotation('rectangle',[0.61, 0.76, 0.03, 0.02],'FaceColor',colorVals(2,:),'EdgeColor','none')

annotation('textbox',[0.64, 0.787, 0.1, 0.1],'String','Slow Gamma','FontSize',20,'EdgeColor','none')
annotation('textbox',[0.64, 0.757, 0.1, 0.1],'String','Fast Gamma','FontSize',20,'EdgeColor','none')
annotation('textbox',[0.64, 0.727, 0.1, 0.1],'String','Slow Gamma Overlap','FontSize',20,'EdgeColor','none')
annotation('textbox',[0.64, 0.697, 0.1, 0.1],'String','Fast Gamma Overlap','FontSize',20,'EdgeColor','none')

annotation('textbox',[0.05,0.95, 0, 0], 'string', 'A','FontSize',30,'FontWeight','bold')