function plotpsdGrid(data,goodElectrodes,fs,ntapers,bandwidth)
%inputs 
% data in form channelsxtimextrials
% goodElectrodes - vector of electrodes with highRMS
% fs - sampling frequency
% ntapers and bandwidth - parameters for chronux
% baseline and event - indices of timepoints with baseline and event - specified by default to 750ms
%%
baseline = 197:1696;
event = 1697:3196;


% create grid layout
gridLayout = flipud(fliplr(reshape(1:81,[9,9])));

dataev = data(:,event,:);
databl = data(:,baseline,:);
%set params for chronux
params.Fs = fs;
% params.err = [2 0.025];
params.tapers = [length(dataev)/params.Fs*bandwidth ntapers];
% params.pad = -1;


for ind = 1:size(dataev,3)
[Sev(:,:,ind),fev] = mtspectrumc(dataev(:,:,ind)',params);
[Sbl(:,:,ind),fbl] = mtspectrumc(databl(:,:,ind)',params);
end
Sev = mean(Sev,3);
Sbl = mean(Sbl,3);
Scorr = 10*log10(Sev./Sbl);

subplotGrid = reshape(1:81,[9,9])';
% plot delta psd's in a grid
for ind = 1:size(Scorr,2)
    elecPos = find(gridLayout==goodElectrodes(ind));
    subplot(9,9,elecPos)
%      plot(fev,Scorr(:,ind),'Color','red','LineWidth',1)
    hold on
    plot(fbl,10*log10(Sbl(:,ind)),'--','Color','black','LineWidth',0.5)
    hold on
    plot(fev,10*log10(Sev(:,ind)),'Color','blue','LineWidth',1)
    title(['Electrode ',num2str(goodElectrodes(ind))])
    xlim([0 150])
end
%%
end