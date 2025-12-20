function plotWavesAllTrials(waveVector,timeVals,boundries,overlap)
if nargin<4
    overlap = [];
end
    colorVals = cat(1,[52 148 186]./255,[236 112  22]./255);
    numTrials = size(waveVector,1);
    trialId = 1:numTrials;
    waveVector(~isnan(waveVector)) = 1;
    plot(timeVals,waveVector(:,:,1).*trialId'-0.3,'LineWidth',1.2,'Color',colorVals(1,:))
    hold on
    if size(waveVector,3)>1
        plot(timeVals,waveVector(:,:,2).*trialId'-0.5,'LineWidth',1.2,'Color',colorVals(2,:))
    end
    yline(trialId,'--','LineWidth',0.5,'Color','black')
    xlim(boundries)
    ylim([1 numTrials])

if ~isempty(overlap)
    plot(timeVals,overlap.*trialId'-0.4,'LineWidth',1,'Color','black')
end

end

    
    
    
    
    
    