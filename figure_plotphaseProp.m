%% load data
[allDataSG,goodElectrodes,timeVals,~,~] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sizePos,orientation,2);
[allDataFG,~,~,~,~] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sizePos,orientation,3);
%% get TW params for both bands
allDataSG = permute(allDataSG,[1,3,2]);
allDataFG = permute(allDataFG,[1,3,2]);

trialNo = 1;
req = 2;
nPerm = 250;

freqs = [25 45];
outputsSG = getTWCircParams(allDataSG(:,:,trialNo),timeVals,goodElectrodes,freqs,req,nPerm);
pgdSG = outputsSG.pgd;

freqs = [45 75];
outputsSG = getTWCircParams(allDataFG(:,:,trialNo),timeVals,goodElectrodes,freqs,req,nPerm);
pgdFG = outputsFG.pgd;
%% time frequency plots of a singel channel, averaged across all trials
chanNo = 16;
movingwin = [0.25 0.025];
Fs = 2000;
params.tapers   = [1 1];
params.pad      = -1;
params.Fs       = Fs;
params.trialave = 0; %averaging across trials
fRange = [0 100];
blRange = [-0.5 0];

%get TF params for slow gamma
data = squeeze(allDataSG(chanNo,:,:));
[S,timeTF,freqVals] = mtspecgramc(data',movingwin,params);
xValToPlot = timeTF+timeVals(1)-1/Fs;
TFPow = log10(abs(S).^2); 
blPos = intersect(find(xValToPlot>=blRange(1)),find(xValToPlot<blRange(2)));
blPower = mean(TFPow(blPos,:,:),1);
logSBL = repmat(blPower,length(xValToPlot),1);
corrPowA(:,:,1) = 10*(mean(TFPow,3)-mean(logSBL,3));

%get TF params for slow gamma
data = squeeze(allDataFG(chanNo,:,:));
[S,timeTF,freqVals] = mtspecgramc(data',movingwin,params);
xValToPlot = timeTF+timeVals(1)-1/Fs;
TFPow = log10(abs(S).^2); 
blPos = intersect(find(xValToPlot>=blRange(1)),find(xValToPlot<blRange(2)));
blPower = mean(TFPow(blPos,:,:),1);
logSBL = repmat(blPower,length(xValToPlot),1);
corrPowA(:,:,2) = 10*(mean(TFPow,3)-mean(logSBL,3));

subplot(2,2,1)
pcolor(xValToPlot,freqVals,corrPowS(:,:,1)');
shading 'interp'
colormap 'jet'
axis([timeRange fRange]);
title(['Time Frequency Map (Slow Gamma) - Electrode',num2str(chanNo)])
hold on
yline(25,'--','Color','black','LineWidth',0.8)
hold on
yline(45,'--','Color','black','LineWidth',0.8)


subplot(2,2,2)
pcolor(xValToPlot,freqVals,corrPowS(:,:,2)');
shading 'interp'
colormap 'jet'
axis([timeRange fRange]);
hold on
yline(45,'--','Color','black','LineWidth',0.8)
hold on
yline(75,'--','Color','black','LineWidth',0.8)
title(['Time Frequency Map (Fast Gamma) - Electrode',num2str(chanNo)])

%% plot the averaged power time series 
% for slow gamma time series
data = mean(abs(hilbert(allDataSG(:,:,trialNo)')),2)';
% color_area = [243 169 114]./255;    % Orange theme
color_line = [236 112  22]./255;
alpha = 0.5;
curve1 = mean(data)+std(data);
curve2 = mean(data)-std(data);
xVals = [timeVals, fliplr(timeVals)];
inBetween = [curve1, fliplr(curve2)];

subplot(2,2,3)
patch = fill(xVals, inBetween, 'g');
set(patch, 'edgecolor', 'none');
set(patch, 'FaceAlpha', alpha);
hold on;
plot(timeVals, mean(data), 'color', color_line,'LineWidth', 2);
hold on
xline(0.26,'--','Color','black','LineWidth',1)
hold on
xline(0.28,'--','Color','black','LineWidth',1)
hold on
xline(0.65,'--','Color','red','LineWidth',1)
hold on
xline(0.675,'--','Color','red','LineWidth',1)
hold on
plot(timeVals,pgd,'*','Color','black')
hold off
title('Slow Gamma Power time series (25-45Hz) with significant PGD')
xlabel('Time(s)')
xlim([timeVals(1),timeVals(end)])


% for fast gamma time series
data = mean(abs(hilbert(allDataFG(:,:,trialNo)')),2)';
% color_area = [128 193 219]./255;    % Blue theme
color_line = [ 52 148 186]./255;
alpha = 0.5;

curve1 = mean(data)+std(data);
curve2 = mean(data)-std(data);
xVals = [timeVals, fliplr(timeVals)];
inBetween = [curve1, fliplr(curve2)];

subplot(2,2,4)
patch = fill(xVals, inBetween, 'g');
set(patch, 'edgecolor', 'none');
set(patch, 'FaceAlpha', alpha);
hold on;
plot(timeVals, mean(data), 'color', color_line,'LineWidth', 2);
hold on
xline(0.26,'--','Color','black','LineWidth',1)
hold on
xline(0.28,'--','Color','black','LineWidth',1)
hold on
xline(0.65,'--','Color','red','LineWidth',1)
hold on
xline(0.675,'--','Color','red','LineWidth',1)
hold on
plot(timeVals,pgd,'*','Color','black')
hold off
title('Fast Gamma Power time series (25-45Hz) with significant PGD')
xlabel('Time(s)')
xlim([timeVals(1),timeVals(end)])
%% plot phase propagation plots - plotted as a seperate figure, combined in photoscape
gridLayout = rot90(reshape(1:81,[9,9]),2); %set the grid layout

data1 = allDataSG(:,:,trialNo);
data2 = allDataFG(:,:,trialNo);

phiGridSG = nan(size(gridLayout,1),size(gridLayout,1),length(data1));
phiGridFG = nan(size(gridLayout,1),size(gridLayout,1),length(data2));

for ind = 1:numel(goodElectrodes)
    [xloc,yloc] = find(gridLayout==goodElectrodes(ind));
    phiGridSG(xloc,yloc,:) = angle(hilbert(data1(ind,:)));
    phiGridFG(xloc,yloc,:) = angle(hilbert(data2(ind,:)));
end
clear data1 data2 
[X,Y] = meshgrid(1:1:9);
timeValuesSG = dsearchn(timeVals',0.26):5:dsearchn(timeVals',0.28);
timeValuesFG = dsearchn(timeVals',0.66):5:dsearchn(timeVals',0.68);

%plot for slow gamma
directionsSG = outputsSG.direction;
for i = 1:length(timeValuesSG)
subplot(3,3,timeValuesSG(i))
imagesc(cos(phiGridFG(:,:,timeValuesSG(i))))
xticks(1)
xticklabels([])
yticks(1)
yticklabels([])
title(['Time:',num2str(timeVals(timeValuesSG(i))),'s'])
caxis([-1 1])
axis square
hold on
U = real(exp(1i*directionsSG(timeValuesSG(i)))); 
V = imag(exp(1i*directionsSG(timeValuesSG(i))));
quiver(X,Y,U,V,'Color','white','LineWidth',1,'AutoScaleFactor',0.5,'parent',hGridPlots(3))                 
hold off
end

%plot for slow gamma
directionsSG = outputsFG.direction;
for i = 1:length(timeValuesSG)
subplot(3,3,timeValuesSG(i))
imagesc(cos(phiGridFG(:,:,timeValuesFG(i))))
xticks(1)
xticklabels([])
yticks(1)
yticklabels([])
title(['Time:',num2str(timeVals(timeValuesFG(i))),'s'])
caxis([-1 1])
axis square
hold on
U = real(exp(1i*directionsFG(timeValuesFG(i)))); 
V = imag(exp(1i*directionsFG(timeValuesFG(i))));
quiver(X,Y,U,V,'Color','white','LineWidth',1,'AutoScaleFactor',0.5,'parent',hGridPlots(3))                 
hold off
end
