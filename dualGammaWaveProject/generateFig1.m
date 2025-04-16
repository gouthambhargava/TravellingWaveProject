%% generate Fig 1 
% load and filter data
trialNo = 1;
elecFrac = 0.5;
selectedElec = 41;
dataPath = 'G:\monkeyData\data';
gridType = 'Microelectrode';
subjectName='alpaH'; expDate = '210817'; protocolName = 'GRF_002'; 
sPos = 2; % spatial frequency
oriPos = 4; % orientation[m1Data,goodElectrodes1,timeVals,rfData1,parameters1] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sPos,oriPos);
freqRangeList{1} = [20 35]; freqRangeList{2} = [40 60];
% get burst data 
% burst calculation parameters
thresholdFactor = 3;
baselinePeriodS = [-0.5 0];
stimulusPeriodS = [0.25 0.75];
analysisPeriodS = [-0.5 1];
filterOrder = 4;
numFrequencyRanges = numel(freqRangeList);

burstTS = zeros(numel(goodElectrodes1),size(m1Data,2),size(m1Data,3),numFrequencyRanges);
filteredSignal = zeros(numel(goodElectrodes1),size(m1Data,2),size(m1Data,3),numFrequencyRanges);
for iFreq=1:numFrequencyRanges
    for iElec=1:numel(goodElectrodes1)
        [~,~,~,burstTS(iElec,:,:,iFreq),filteredSignal(iElec,:,:,iFreq)] = getHilbertBurst(squeeze(m1Data(iElec,:,:)),timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,freqRangeList{iFreq},filterOrder,1,analysisPeriodS);
     end
end

% generate figure
rfColors = jet(numel(goodElectrodes1));
figure(2)
tiledlayout(6,2,'TileSpacing','Compact','Padding','Loose')
grid1 = nexttile(1,[2,1]);
grid2 = nexttile(2,[2,1]);
grid3 = nexttile(5);
grid4 = nexttile(6);
grid5 = nexttile(7,[2,1]);
grid6 = nexttile(8,[2,1]);
grid7 = nexttile(11);
grid8 = nexttile(12);

colorNames = gray(3);
axisLims = {[1 80],[-0.5 1],[-2 2]};
plotTFMT(grid1,m1Data(selectedElec,:,:),timeVals,axisLims,freqRangeList,colorNames,'delta',[],[0 0.8])
plotTFMT(grid2,m1Data(selectedElec,:,:),timeVals,axisLims,freqRangeList,colorNames,'single',trialNo,[0 0.8])
clim(grid1,[-3 5])
ylabel(grid1,'Frequency')
c = colorbar(grid1,'northoutside');
colorbar(grid2,'northoutside')

% plot time series with bursts
colorVals = cat(1,[52 148 186]./255,[236 112  22]./255);
plot(grid3,timeVals,squeeze(filteredSignal(selectedElec,trialNo,:,1)),'LineWidth',1,'Color',colorVals(1,:));
hold(grid3,'on')
burstInt = squeeze(burstTS(selectedElec,trialNo,:,1));
burstInt(burstInt==1) = 0;
boundsSG = getLims(burstInt,1);
plot(grid3,timeVals,burstInt,'LineWidth',1,'Color','black');
xlim(grid3,[-0.5 1])
ylim(grid3,[-50 50])
% xlabel(grid3,'Time(s)')
ylabel(grid3,'Voltage(\muV)')
title(grid3,'Slow Gamma(20-35Hz)')
xline(grid3,[0 0.8],'LineWidth',1)

plot(grid4,timeVals,squeeze(filteredSignal(selectedElec,trialNo,:,2)),'LineWidth',1,'Color',colorVals(2,:));
hold(grid4,'on')
burstInt = squeeze(burstTS(selectedElec,trialNo,:,2));
burstInt(burstInt==1) = 0;
boundsFG = getLims(burstInt,1);
plot(grid4,timeVals,burstInt,'LineWidth',1,'Color','black');
xlim(grid4,[-0.5 1])
ylim(grid4,[-50 50])
% ylabel(grid4,'Voltage(\muV)')
% xlabel(grid4,'Time(s)')
title(grid4,'Fast Gamma(40-60Hz)')
xline(grid4,[0 0.8],'LineWidth',1)

% plot burst locs 
chans1 = 1:numel(goodElectrodes1);
plot(timeVals,squeeze(burstTS(:,trialNo,:,1)).*chans1','.','color',colorVals(1,:),'linewidth',1,'parent',grid5);
xlim(grid5,[-0.5 1])
ylabel(grid5,'Electrodes')
% xlabel(grid5,'Time(s)') 
hold(grid5,'on')
ylim(grid5,[0 77])
plot(timeVals,squeeze(burstTS(selectedElec,trialNo,:,1)).*chans1(selectedElec)','.','color',rfColors(1,:),'linewidth',1,'parent',grid5);
xline(grid5,[0 0.8],'LineWidth',1)

chans1 = 1:numel(goodElectrodes1);
plot(timeVals,squeeze(burstTS(:,trialNo,:,2)).*chans1','.','color',colorVals(2,:),'linewidth',1,'parent',grid6);
xlim(grid6,[-0.5 1])
% ylabel(grid6,'Electrodes')
% xlabel(grid6,'Time(s)') 
hold(grid6,'on')
ylim(grid6,[0 77])
plot(timeVals,squeeze(burstTS(selectedElec,trialNo,:,2)).*chans1(selectedElec)','.','color',rfColors(1,:),'linewidth',1,'parent',grid6);
xline(grid6,[0 0.8],'LineWidth',1)

% plot the burst fraction data
burstFrac1 = squeeze(nansum(squeeze(burstTS(:,trialNo,:,:)),1));
ylineLim1 = elecFrac*numel(goodElectrodes1);
fractionPlot1 = burstFrac1;
fractionPlot1(fractionPlot1<77*elecFrac) = nan;
fractionPlot1(fractionPlot1>0) = 0;

plot(grid7,timeVals,burstFrac1(:,1),'LineWidth',1.5,'Color',colorVals(1,:))
hold(grid7,'on')
% yline(ylineLim1,'--','LineWidth',1.5,'Color','k')
plot(grid7,timeVals,ones(1,numel(timeVals))*ylineLim1,'--','LineWidth',1.5,'Color','k')
plot(grid7,timeVals,fractionPlot1(:,1),'|','Color',colorVals(1,:))
ylim(grid7,[-10 80])
xlabel(grid7,'Time(s)')
xlim(grid7,[-0.5 1])
yticks(grid7,[0 46.2 77])
yticklabels(grid7,{'0','0.5','1'})
ylabel(grid7,'Elec Fraction')
xline(grid7,[0 0.8],'LineWidth',1)

plot(grid8,timeVals,burstFrac1(:,2),'LineWidth',1.5,'Color',colorVals(2,:))
hold(grid8,'on')
yline(ylineLim1,'--','LineWidth',1.5,'Color','k')
plot(grid8,timeVals,fractionPlot1(:,2),'|','Color',colorVals(2,:))
ylim(grid8,[-10 80])
xlabel(grid8,'Time(s)')
xlim(grid8,[-0.5 1])
yticks(grid8,[0 46.2 77])
yticklabels(grid8,{'0','0.5','1'})
% ylabel(grid8,'EF')
xline(grid8,[0 0.8],'LineWidth',1)

% plot lines on the tf plot
hold(grid2,'on')
freqs = [24 24 30 32 24];
for i = 1:size(boundsSG,2)
    plotInt = nan(1,numel(timeVals));
    plotInt(boundsSG(1,i):boundsSG(2,i)) = freqs(i);
    plot(grid2,timeVals,plotInt,'LineWidth',1,'Color','black');
end

freqs = [45 44 44 44 52 48 56];
for i = 1:size(boundsFG,2)
    plotInt = nan(1,numel(timeVals));
    plotInt(boundsFG(1,i):boundsFG(2,i)) = freqs(i);
    plot(grid2,timeVals,plotInt,'LineWidth',1,'Color','black');
end
%
annotation('textbox',[0.02,0.95, 0, 0], 'string', 'A','FontSize',24,'FontWeight','bold')
annotation('textbox',[0.48,0.95, 0, 0], 'string', 'B','FontSize',24,'FontWeight','bold')
annotation('textbox',[0.02,0.65, 0, 0], 'string', 'C','FontSize',24,'FontWeight','bold')
annotation('textbox',[0.48,0.65, 0, 0], 'string', 'D','FontSize',24,'FontWeight','bold')
annotation('textbox',[0.02,0.48, 0, 0], 'string', 'E','FontSize',24,'FontWeight','bold')
annotation('textbox',[0.48,0.48, 0, 0], 'string', 'F','FontSize',24,'FontWeight','bold')
annotation('textbox',[0.02,0.21, 0, 0], 'string', 'G','FontSize',24,'FontWeight','bold')
annotation('textbox',[0.48,0.21, 0, 0], 'string', 'H','FontSize',24,'FontWeight','bold')

