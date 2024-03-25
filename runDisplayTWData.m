% This code runs the main display program (displayTWData)
%% select data of particular orientation
dataPath = 'G:\alpa\data\';
gridType = 'Microelectrode';

subjectName='alpaH'; expDate = '210817'; protocolName = 'GRF_002'; selectedElectrodes = [1 41 81]; freqRangeList{1} = [6 12]; freqRangeList{2} = [45 70];

sPos = 2; % spatial frequency: 0.5 (1), 1(2), 2 (3), 4 (4), 8 (5), all SFs (6). Note that the same code can be used for the size project also later where stimulus size is changed instead of spatial frequency
oriPos = 3; % orientation: 0 (1), 22.5 (2), 45 (3), 67.5 (4), 90 (5), 112.5 (6), 135 (7), 157.5 (8), all orientations (9)

analysisMethod = 'hilbert';
displayTWData(subjectName,expDate,protocolName,dataPath,sPos,oriPos,selectedElectrodes,analysisMethod,freqRangeList);
% 
% %% To get travelling wave params for a single trial
% trialNo = 1;
% freqs = [30 60];
% nPerm = []; %performance similar with 250/500/1000 permutations. Decreases from 250 downwards.
% [data,goodElectrodes,timeVals,rfData,parameters] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sizePos,orientation,4);
% data = squeeze(data(:,trialNo,:));
% TW_outputs = getTWCircParams(data,timeVals,goodElectrodes,freqs,0,nPerm);
% 
% %% plot TW outputs
% pgd = TW_outputs.pgd(:,1);
% pgd(isnan(pgd)) = 0;
% pgdP = TW_outputs.pgdPerm;
% pgdP(isnan(pgdP)) = 0;
% sigValues = find(pgd>pgdP);
% speed = TW_outputs.speed;
% direction = TW_outputs.direction;
% clusterSize = TW_outputs.clusters;
% ampTS = abs(hilbert(data')).^2';
% %%
% subplot(2,3,1:3)
% pgdSig = nan(size(pgd));
% pgdSig(sigValues) = 0;
% plot(timeVals,mean(ampTS),'LineWidth',1,'Color','blue')
% hold on
% plot(timeVals,pgdSig,'--','LineWidth',2,'Color','black')
% hold off
% title('Average Amplitude TS and significant PGD time points')
% xlabel('Time(s)')
% 
% subplot(2,3,4)
% polarhistogram(direction,36,'normalization','pdf','FaceColor',[0.3010 0.7450 0.9330])
% hold on
% directionSig = nan(size(direction));
% directionSig(sigValues) = direction(sigValues);
% polarhistogram(directionSig,36,'normalization','pdf','FaceColor',[0.8500 0.3250 0.0980])
% hold off
% legend('All PGD','Significant PGD')
% title('Direction of propagation for single trial')
% 
% subplot(2,3,5)
% histogram(speed,'normalization','pdf','FaceColor',[0.3010 0.7450 0.9330])
% hold on
% speedSig = nan(size(speed));
% speedSig(sigValues) = speed(sigValues);
% histogram(speedSig,'normalization','pdf','FaceColor',[0.8500 0.3250 0.0980])
% hold off
% xlabel('Speed(m/s')
% legend('All PGD','Significant PGD')
% title('Speed of propagation for single trial')
% 
% subplot(2,3,6)
% histogram(clusterSize,'normalization','pdf','FaceColor',[0.3010 0.7450 0.9330])
% hold on
% clustersSig = nan(size(clusterSize));
% clustersSig(sigValues) = clusterSize(sigValues);
% histogram(clustersSig,'normalization','pdf','FaceColor',[0.8500 0.3250 0.0980])
% hold off
% xlabel('No of electrodes')
% legend('All PGD','Significant PGD')
% title('Histgram of No of electrodes recruited in a single trail')
% 
