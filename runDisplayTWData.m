% This code runs the main display program (displayTWData)

% The following programs should be downloaded and added to Matlab path 
% 1. Chronux: download from http://chronux.org/
% 2. ComonPrograms: https://github.com/supratimray/CommonPrograms
% 3. Gamma Length Project repository: https://github.com/supratimray/GammaLengthProjectCodes)
% 4. Supporting matlab files for MP: http://www.fuw.edu.pl/~durka/software/mp/.
% 5. Cricular statistics toolbox: https://github.com/circstat/circstat-matlab

%% select data of particular orientation
dataPath = 'G:\Bhargava_TravelingWaveProject\data\';
gridType = 'Microelectrode';
subjectName='alpaH'; expDate = '130418'; protocolName = 'GRF_004'; 
sizePos = 7; %0.1500 (1)   0.3000 (2)    0.6000 (3)   1.2000 (4)   2.4000 (5)    4.7900 (6)   9.6000 (7) all (8)
% only 6 sizes (2-7) for kesari and protocols other than the one being used
% for alpaH, adjust accordingly
req = 1; %load reconstructed MP data for whole gamma band
selectedElectrodes = [40 41 42 43 44 45];
orientation = 9; % 0 (1)   22 (2)   45 (3)   67 (4)   90 (5)  112 (6)  135 (7)  157 (8) all orientations (9)
displayTWData(subjectName,expDate,protocolName,dataPath,sizePos,orientation,selectedElectrodes)

% To get travelling wave params for a single trial
trialNo = 1;
nPerm = 250; %performance similar with 250/500/1000 permutations. Decreases from 250 downwards.
[data,goodElectrodes,timeVals,rfData,parameters] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sizePos,orientation,4);
data = squeeze(data(:,trialNo,:));
TW_outputs = getTWCircParams(data,timeVals,goodElectrodes,freqs,0,nPerm);

% plot TW outputs
pgd = TW_outputs.pgd;
pgdP = TW_outputs.pgdPerm;
asigValues = find(pgdP>pgd);
speed = TW_outputs.speed;
direction = TW_outputs.direction;
clusterSize = TW_outputs.clusters;
ampTS = abs(hilbert(data')).^2';

subplot(2,3,1:3)
pgdSig = pgd;
pgdSig(asigValues) = nan;
pgdSig(pgdSig=0) = nan;
pgdSig(isnan(pgdSig)) = 0;
plot(timeVals,mean(ampTS),'LineWidth',1,'Color','blue')
hold on
plot(timeVals,pgdSig,'--','LineWidth',2,'Color','black')
hold off
title('Average Amplitude TS and significant PGD time points')
xlabel('Time(s)')

subplot(2,3,4)
polarhistogram(direction,36,'normalization','pdf','FaceColor',[0.3010 0.7450 0.9330])
hold on
directionSig = direction;
directionSig(asigValues) = nan;
polarhistogram(directionSig,36,'normalization','pdf','FaceColor',[0.8500 0.3250 0.0980])
hold off
legend('All PGD','Significant PGD')
title('Direction of propagation for single trial')

subplot(2,3,5)
histogram(speed,'normalization','pdf','FaceColor',[0.3010 0.7450 0.9330])
hold on
speedSig = speed;
speedSig(asigValues) = nan;
polarhistogram(speedSig,'normalization','pdf','FaceColor',[0.8500 0.3250 0.0980])
hold off
xlabel('Speed(m/s')
legend('All PGD','Significant PGD')
title('Speed of propagation for single trial')

subplot(2,3,6)
histogram(clusters,'normalization','pdf','FaceColor',[0.3010 0.7450 0.9330])
hold on
clustersSig = clusters;
clustersSig(asigValues) = nan;
histogram(clustersSig,'normalization','pdf','FaceColor',[0.8500 0.3250 0.0980])
hold off
xlabel('No of electrodes')
legend('All PGD','Significant PGD')
title('Histgram of No of electrodes recruited in a single trail')

