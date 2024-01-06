function getTwparamsPlot(data,timeVals,goodElectrodes,times,outputs,sigValue,trial)
%% input params
% data - filtered data
% times - time points (x4) of phigrid to be plotted
% output structure obtained from getTWparams which contains: pgd, speed, direction and/or cluster
% sigValue is the threshold of the pgd beyond which it is considered
% significant, by default to be considered 0.5.
%% TW plots
%reshape filtData to grid
gridLayout = flipud(fliplr(reshape(1:81,[9,9]))); %set the grid layout for alpaH
for i = 1:numel(goodElectrodes)
        [x,y] = find(gridLayout==goodElectrodes(i));
        phiGrid(x,y,:) = unwrap(angle(hilbert(data(i,:,trial))));
end 

%get timevals

subplot(4,4,1)
%plot imagesc of the tw propagation 1
freq = '\gamma';
timeStamp = dsearchn(timeVals',times(1));
imagesc(cos(phiGrid(:,:,timeStamp)))
title(['Traveling Wave',freq,'-',num2str(timeVals(timeStamp))])
caxis([min(cos(phiGrid),[],'all'),max(cos(phiGrid),[],'all')])

subplot(4,4,2)
%plot imagesc of the tw propagation 2
timeStamp = dsearchn(timeVals',times(2));
imagesc(cos(phiGrid(:,:,timeStamp)))
caxis([min(cos(phiGrid),[],'all'),max(cos(phiGrid),[],'all')])
title(['Traveling Wave',freq,'-',num2str(timeVals(timeStamp))])

subplot(4,4,3)
%plot imagesc of the tw propagation 3
timeStamp = dsearchn(timeVals',times(3));
imagesc(cos(phiGrid(:,:,timeStamp)))
caxis([min(cos(phiGrid),[],'all'),max(cos(phiGrid),[],'all')])
title(['Traveling Wave',freq,'-',num2str(timeVals(timeStamp))])

subplot(4,4,4)
%plot imagesc of the tw propagation 3
timeStamp = dsearchn(timeVals',times(4));
imagesc(cos(phiGrid(:,:,timeStamp)))
caxis([min(cos(phiGrid),[],'all'),max(cos(phiGrid),[],'all')])
title(['Traveling Wave',freq,'-',num2str(timeVals(timeStamp))])


subplot(4,4,5:8)
%plot the stair plot of the pgd/cluster
if outputs.analysis=='Gradient'
stairs(timeVals,outputs.pgd(:,trial),'LineWidth',1,'Color','red')
xlim([timeVals(1),timeVals(end)])
title(['Phase Gradient Directionality for trial-',num2str(trial)])
else
stairs(timeVals,outputs.cluster(:,trial),'LineWidth',1,'Color','red')
xlim([timeVals(1),timeVals(end)])
title(['Cluster size of travelling wave for trial-',num2str(trial)])
end
hold on
sigPGDval = find(outputs.pgd>sigValue); %get significant pgd values
if ~isempty(sigPGDval)
ypts = 0.7*ones(1,length(sigPGDval));
plot(timeVals(values),ypts,'.')
ylim([0 0.8])
yticks(sigValue)
yticklabels({'Significant'})
end
xlabel('Time(s)')
hold off

subplot(4,4,9:12)
%plot the average power time series with sig PGD values
powData = abs(hilbert(mean(data(:,:,trial)))).^2;
stairs(timeVals,powData,'LineWidth',1,'Color','blue')
xlim([timeVals(1),timeVals(end)])
hold on
sigPGDval = find(outputs.pgd>sigValue); %get significant pgd values
if ~isempty(sigPGDval)
ypts = 0.7*ones(1,length(sigPGDval));
plot(timeVals(values),ypts,'.')
ylim([0 0.8])
yticks(sigValue)
yticklabels({'Significant'})
end
xlabel('Time(s)')
hold off
title(['Phase Gradient Directionality for trial-',num2str(trial)])



subplot(4,4,14)
%plot polar histogram of the direction
polarhistogram(outputs.direction(:,trial))
hold on
if ~isempty(sigPGDval)
polarhistogram(direction(sigPDGVal,trial))
labels('All times','Only Significant')
end
title('Direction distribution')


subplot(4,4,15)
%plot polar histograms of direction seperately for baseline and stimulus
histogram(outputs.speed(:,trial),'Normalization','pdf')
hold on
if ~isempty(sigPGDval)
histogram(speed(sigPGDval,trial),'Normalization','pdf')
hold on
legend('All Timepoints','Significant')
end
title('Traveling wave speed')
xlabel('Speed (m/s)')
ylabel('PDF')
end



