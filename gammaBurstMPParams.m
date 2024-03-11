function [gammaAtom, sigBursts, burstTS, gaborInfo] = gammaBurstMPParams(lengthList,timeList,freqList,gaborInfo,timeVals)
%Extract params from MP gamma burst algo for TW calculation and plotting.
%exclude bursts that lasted over 0.8s and select the burst with the max
%length from the remaining
gammaAtom = [];
sigBursts = cat(2,lengthList,timeList,freqList);
sigBursts(sigBursts(:,1)>0.8,:) = [];
sigBursts = sigBursts(sigBursts(:,1)==max(sigBursts(:,1)),:);

burstTS = zeros(1,length(timeVals));
burstTS(intersect(find(timeVals>(sigBursts(2)-sigBursts(1)/2)),find(timeVals<(sigBursts(2)+sigBursts(1)/2)))) = 1;

gammaAtom(chanNo) = find(squeeze(gaborInfo(trialNo,:,2))==burstTS(3));
gaborInfo = gaborInfo(gammaAtom,:);
end
