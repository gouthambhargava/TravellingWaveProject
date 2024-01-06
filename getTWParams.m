function TWimagePalette(dataPath,subjectName,gridType,data,goodElectrodes,timeVals,timeFrame,videoTitle,trial,stimSize) 
%% inputs
% dataPath = 'G:\Bhargava_TravelingWaveProject\data\'; 
% subjectName='alpaH'; gridType = 'Microelectrode';
% timeFrame - vector of values indicating the segment of data to be
% visualized, can be specified in seconds as [0 1] which according to timeVals is 1696
% and 3696
timeFrame = dsearchn(timeVals',timeFrame(1)):dsearchn(timeVals',timeFrame(2));
writerObj = VideoWriter(videoTitle); %open the video writer object

gridLayout = flipud(fliplr(reshape(1:81,[9,9]))); %set grid layout alpaH
ampGrid = nan(size(gridLayout,1),size(gridLayout,2),size(data,2));
phiGrid = nan(size(gridLayout,1),size(gridLayout,2),size(data,2));

for i = 1:numel(goodElectrodes)
    [x,y] = find(gridLayout==goodElectrodes(i));
    ampGrid(x,y,:) = abs(hilbert(data(i,:,trial))).^2; %power grid 
    phiGrid(x,y,:) = unwrap(angle(hilbert(data(i,:,trial)))); %phase grid
end
clear x y

[goodDPos,distanceBins] = getStimDistBins(dataPath,subjectName,gridType,stimSize);% bin the electrodes according to distance
colorVals = hsv(numel(distanceBins));
%average electrodes in each bin
distPlotTS = zeros(numel(goodDPos),size(data,2));
for i = 1:numel(goodDPos)
    distPlotTS(i,:) = mean(data(goodDPos{i},:,trial));
end

frames = [];
vals = timeVals(timeFrame);
figure('units','normalized','outerposition',[0 0 1 1])
for ind = 1:length(vals)
    
    subplot(6,4,[1,2,5,6])
    %plot the stimulus size and location
    plotStimLoc(dataPath,subjectName,gridType,stimSize)

    subplot(6,4,[9,10,13,14])
    %plot power in the grid with superimposed phase 
    data1 = ampGrid(:,:,timeFrame);
    data2 = phiGrid(:,:,timeFrame);
    [X,Y] = meshgrid(1:1:9);
    U = cos(data2(:,:,ind)); 
    V = sin(data2(:,:,ind)); 
    imagesc(data1(:,:,ind))
    xticks([2 4 6 8])
    xticklabels({'','','',''})
    yticks([1 2 3 4 5 6 7 8 9])
    yticklabels({'','','','','','','','',''})
    caxis([min(data1,[],'all'),max(data1,[],'all')])
    hold on
    quiver(X,Y,U,V,'Color','white','LineWidth',1,'AutoScaleFactor',0.5)
    hold off
    title(['Power across grid at time:',num2str(vals(ind)),'s'])
    axis square
    colormap(hot)
    
    subplot(6,4,[17,18,21,22])
    data1 = phiGrid(:,:,timeFrame);
    imagesc(cos(data1(:,:,ind)))
    xticks([2 4 6 8])
    xticklabels({'','','',''})
    yticks([1 2 3 4 5 6 7 8 9])
    yticklabels({'','','','','','','','',''})
    caxis([min(cos(data1),[],'all'),max(cos(data1),[],'all')])
    title(['Phase propagation at time:',num2str(vals(ind)),'s'])
    axis square  
    colormap(hot)
    
    subplot(6,4,3:4)
    plot(vals,abs(hilbert(mean(distPlotTS(:,timeFrame)))).^2,'LineWidth',1.5,'Color','black')
    hold on
    xline(vals(ind))
    hold off
    title(['Gamma power time series:',num2str(vals(ind))])
    
    subplot(6,4,7:8)
    plot(vals,distPlotTS(1,timeFrame),'LineWidth',1.5,'Color',colorVals(1,:))
    hold on
    xline(vals(ind))
    hold off
    title(['Averaged filtered LFP at distance:',num2str(distanceBins(1))])
    
    subplot(6,4,11:12)
    plot(vals,distPlotTS(2,timeFrame),'LineWidth',2,'Color',colorVals(2,:))    
    hold on
    xline(vals(ind))
    hold off
    title(['Averaged filtered LFP at distance:',num2str(distanceBins(2))])
    
    subplot(6,4,15:16)
    plot(vals,distPlotTS(3,timeFrame),'LineWidth',2,'Color',colorVals(3,:))    
    hold on
    xline(vals(ind))
    hold off
    title(['Averaged filtered LFP at distance:',num2str(distanceBins(3))])
    
    subplot(6,4,19:20)
    plot(vals,distPlotTS(4,timeFrame),'LineWidth',2,'Color',colorVals(4,:))    
    hold on
    xline(vals(ind))
    hold off
    title(['Averaged filtered LFP at distance:',num2str(distanceBins(4))])
    
    subplot(6,4,23:24)
    plot(vals,distPlotTS(5,timeFrame),'LineWidth',2,'Color',colorVals(5,:))    
    hold on
    xline(vals(ind))
    hold off
    title(['Averaged filtered LFP at distance:',num2str(distanceBins(5))])
    frames = [frames getframe(gcf)];
end
%%


writerObj.FrameRate = 30; % set the seconds per image
% open the video writer
open(writerObj);
for i=1:length(frames)
    % convert the image to a frame
    frame = frames(i) ;    
    writeVideo(writerObj, frame);
end
close(writerObj);
end
