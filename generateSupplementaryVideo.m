%% get video TW
%% for alpa
sPos = 2;
oriPos = 4;
req = 2;
dataPath = 'F:\monkeyData\data';
gridType = 'Microelectrode';
freqRangeList{1} = [20 35]; freqRangeList{2} = [40 60];
if req ==1
    subjectName='alpaH'; expDate = '210817'; protocolName = 'GRF_002'; 
    [mData,goodElectrodes,timeVals,~,~] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sPos,oriPos);
outputsTW = outputsTWA24;
else
    % load data for M2
    subjectName='kesariH'; expDate = '270218'; protocolName = 'GRF_001';  
    [mData,goodElectrodes,timeVals,~,~] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sPos,oriPos);
outputsTW = outputsTWK23;
end

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
trial = 15;
wave = 1;
wobbleReq = 2;
[newBounds,sigCells,allDirSg,allDirFg,allUniqueDirs] = getOverlappingWaves(outputsTW,timeVals,0.5,wobbleReq);
[~,~,~,~,~,~,pgd1] = getWaveSegments(outputsTW(1,:),timeVals,wobbleReq);
[~,~,~,~,~,~,pgd2] = getWaveSegments(outputsTW(2,:),timeVals,wobbleReq);


pgd1 = pgd1(sigCells,:);
pgd2 = pgd2(sigCells,:);
duration1 = newBounds{1,trial}(:,wave);
duration2 = newBounds{2,trial}(:,wave);
durOverlap = intersect(duration1(1):duration1(2),duration2(1):duration2(2));
% durOverlap for alpa for trail 15/1 has been taken as 2708:2744
durOverlap = 2708:2775;

gridLayout = rot90(reshape(1:81,[9,9]),2); %set the grid layout
phaseData1 = nan(size(gridLayout,1),size(gridLayout,2),numel(timeVals));
phaseData2 = nan(size(gridLayout,1),size(gridLayout,2),numel(timeVals));
burstFrac = squeeze(nansum(squeeze(burstTS(:,trial,:,:)),1));
burstFrac(burstFrac<0.5*numGoodElectrodes) = nan;
burstFrac = burstFrac/numGoodElectrodes;
pgd1 = pgd1(trial,:);
pgd2 = pgd2(trial,:);
for i = 1:numel(goodElectrodes)
    [x,y] = find(gridLayout==goodElectrodes(i));
    phaseData1(x,y,:) = angle(hilbert(squeeze(filteredSignal(i,trial,:,1))));
    phaseData2(x,y,:) = angle(hilbert(squeeze(filteredSignal(i,trial,:,2))));
end


phaseData1 = phaseData1(:,:,durOverlap(1):durOverlap(end));
phaseData2 = phaseData2(:,:,durOverlap(1):durOverlap(end));


timeValues = timeVals(durOverlap(1):durOverlap(end));
u1 = cos(allDirSg(trial,durOverlap(1):durOverlap(end)));
v1 = sin(allDirSg(trial,durOverlap(1):durOverlap(end)));
u2 = cos(allDirFg(trial,durOverlap(1):durOverlap(end)));
v2 = sin(allDirFg(trial,durOverlap(1):durOverlap(end)));
%
f = figure;
set(gcf, 'Position', get(0, 'Screensize'));
% f.Position = [100 100 540 400];
colorVals = cat(1,[52 148 186]./255,[236 112  22]./255);
[x, y] = meshgrid(1:9,1:9);
pgd1(isnan(pgd1)) = 0;
pgd2(isnan(pgd2)) = 0;

a = allDirSg(trial,:);
b = allDirFg(trial,:);
c = nan(1,numel(timeVals));
c(durOverlap) = 1;
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
        writerObj = VideoWriter('alpaTest2.avi');
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
    quiver(x,y,v1(i).*ones(size(x)),u1(i).*ones(size(y)),'color','white','AutoScaleFactor',0.9,'parent',phaseGrid1)
    axis(phaseGrid1,'square')
    axis(phaseGrid1, 'off')
    title(phaseGrid1,['Slow Gamma:', num2str(timeValues(i)),'s'])
    colormap(phaseGrid1,'jet')
    imagesc(cos(phaseData2(:,:,i)),'Parent',phaseGrid2)
    clim(phaseGrid2,[-1 1])
    hold(phaseGrid2,'on')
    quiver(x,y,v2(i).*ones(size(x)),u2(i).*ones(size(y)),'color','white','AutoScaleFactor',0.9,'parent',phaseGrid2)
    axis(phaseGrid2,'square')
    axis(phaseGrid2, 'off')
    title(phaseGrid2,['Fast Gamma:', num2str(timeValues(i)),'s'])
    colormap(phaseGrid2,'jet')
    pause(0.1)
    frames = getframe(gcf);
    writeVideo(writerObj, frames);
    delete(h)
end
close(writerObj)