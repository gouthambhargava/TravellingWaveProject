function [burstTS, filteredSignal] = runHilbertBurstLength(data,timeVals,freqs,req)

filterOrder = 4; 
thresholdFactor = 1;
baselinePeriodS = [-0.5 0]; 
stimulusPeriodS = [0.25 0.75];

    for i = 1:size(data,1)
        [~,~,~,burstTS(i,:),filteredSignal(i,:),~] = getBurstLengthHilbert(data(i,:),timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,freqs,filterOrder,req);
    end
    burstTS(isnan(burstTS)) = 0;
end