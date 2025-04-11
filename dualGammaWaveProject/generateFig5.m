%% generate fig 5
req = 1;
minBurstSize = 50;
wobble = 10;
sPos = 2;
oriPos = 1:8;
thresh = 0.5;
binEdges = 0:0.05:1;
outputsTW1 = {outputsTWA21,outputsTWA22,outputsTWA23,outputsTWA24,outputsTWA25,outputsTWA26,outputsTWA27,outputsTWA28};
outputsTW2 = {outputsTWK21,outputsTWK22,outputsTWK23,outputsTWK24,outputsTWK25,outputsTWK26,outputsTWK27,outputsTWK28};

for i = 1:numel(oriPos)
[sgBursts1{i},fgBursts1{i}] = getBurstOverlap(outputsTW1{i},sPos,i,thresh,minBurstSize,binEdges,req);
end
%
req = 2;
for i = 1:numel(oriPos)
[sgBursts2{i},fgBursts2{i}] = getBurstOverlap(outputsTW2{i},sPos,i,thresh,minBurstSize,binEdges,req);
end

%%
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




