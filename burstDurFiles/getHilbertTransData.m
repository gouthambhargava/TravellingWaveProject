function [ampGrid,phiGrid,mag] = getHilbertTransData(signal,timeVals,gammaFreqRangeHz,filtOrder,goodElectrodes)
%signal is a m by n matrix where m is the time points n is the number of
%channels
%check dimensions of signal and reshape if required
signal = squeeze(signal);
if size(signal,1)<size(signal,2)
    signal = signal';
end
Fs=1/(timeVals(2)-timeVals(1));
normBand=gammaFreqRangeHz/(Fs/2);
[b,a]=butter(filtOrder,normBand,'bandpass');
bpfSignal=filtfilt(b,a,signal);

% phiVal = angle(hilbert(bpfSignal))';
% powVal = abs(hilbert(bpfSignal)).^2';

% get phi and amp grid
gridLayout = flipud(fliplr(reshape(1:81,[9,9]))); %set grid layout alpaH
ampGrid = nan(size(gridLayout,1),size(gridLayout,2),size(bpfSignal,1));
phiGrid = nan(size(gridLayout,1),size(gridLayout,2),size(bpfSignal,1));

   for i = 1:numel(goodElectrodes)
       [x,y] = find(gridLayout==goodElectrodes(i));
       ampGrid(x,y,:) = abs(hilbert(bpfSignal(:,i))).^2;
       phiGrid(x,y,:) = angle(hilbert(bpfSignal(:,i)));
   end
  
   %calculate phase coherence 
   %add the phases together and calculate the magnitude of this for each
   %time point.
   mag = zeros(1,length(phiGrid));
   for ind = 1:size(phiGrid,3)
       phiValues = reshape(phiGrid(:,:,ind),[1,numel(phiGrid(:,:,ind))]);
       phiValues(isnan(phiValues)) = [];
       mag(ind) = sqrt(sum(phiValues).^2);     
   end
%    ampGrid = ampGrid(:,:,timeRange);
%    phiGrid = phiGrid(:,:,timeRange);
   clear x y
   
end