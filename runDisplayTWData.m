% This code runs the main display program (displayTWData)
%% select data of particular orientation
dataPath = 'E:\IISc_exp\IISC_work\GitScripts\monkeyData\data\';
gridType = 'Microelectrode';
subjectName='alpaH'; expDate = '130418'; protocolName = 'GRF_004'; 
sizePos = 7; %0.1500 (1)   0.3000 (2)    0.6000 (3)   1.2000 (4)   2.4000 (5)    4.7900 (6)   9.6000 (7) all (8)
% only 6 sizes (2-7) for kesari and protocols other than the one being used
% for alpaH, adjust accordingly
req = 4; %load reconstructed MP data for whole gamma band
selectedElectrodes = [40 41 42 43 44 45];
orientation = 9; % 0 (1)   22 (2)   45 (3)   67 (4)   90 (5)  112 (6)  135 (7)  157 (8) all orientations (9)
displayTWData(subjectName,expDate,protocolName,dataPath,sizePos,orientation,selectedElectrodes,req)

%% To get travelling wave params for a single trial
trialNo = 1;
freqs = [30 60];
nPerm = []; %performance similar with 250/500/1000 permutations. Decreases from 250 downwards.
[data,goodElectrodes,timeVals,rfData,parameters] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sizePos,orientation,4);
data = squeeze(data(:,trialNo,:));
TW_outputs = getTWCircParams(data,timeVals,goodElectrodes,freqs,0,nPerm);

%% plot TW outputs
pgd = TW_outputs.pgd(:,1);
pgd(isnan(pgd)) = 0;
pgdP = TW_outputs.pgdPerm;
pgdP(isnan(pgdP)) = 0;
sigValues = find(pgd>pgdP);
speed = TW_outputs.speed;
direction = TW_outputs.direction;
clusterSize = TW_outputs.clusters;
ampTS = abs(hilbert(data')).^2';
%%
subplot(2,3,1:3)
pgdSig = nan(size(pgd));
pgdSig(sigValues) = 0;
plot(timeVals,mean(ampTS),'LineWidth',1,'Color','blue')
hold on
plot(timeVals,pgdSig,'--','LineWidth',2,'Color','black')
hold off
title('Average Amplitude TS and significant PGD time points')
xlabel('Time(s)')

subplot(2,3,4)
polarhistogram(direction,36,'normalization','pdf','FaceColor',[0.3010 0.7450 0.9330])
hold on
directionSig = nan(size(direction));
directionSig(sigValues) = direction(sigValues);
polarhistogram(directionSig,36,'normalization','pdf','FaceColor',[0.8500 0.3250 0.0980])
hold off
legend('All PGD','Significant PGD')
title('Direction of propagation for single trial')

subplot(2,3,5)
histogram(speed,'normalization','pdf','FaceColor',[0.3010 0.7450 0.9330])
hold on
speedSig = nan(size(speed));
speedSig(sigValues) = speed(sigValues);
histogram(speedSig,'normalization','pdf','FaceColor',[0.8500 0.3250 0.0980])
hold off
xlabel('Speed(m/s')
legend('All PGD','Significant PGD')
title('Speed of propagation for single trial')

subplot(2,3,6)
histogram(clusterSize,'normalization','pdf','FaceColor',[0.3010 0.7450 0.9330])
hold on
clustersSig = nan(size(clusterSize));
clustersSig(sigValues) = clusterSize(sigValues);
histogram(clustersSig,'normalization','pdf','FaceColor',[0.8500 0.3250 0.0980])
hold off
xlabel('No of electrodes')
legend('All PGD','Significant PGD')
title('Histgram of No of electrodes recruited in a single trail')

