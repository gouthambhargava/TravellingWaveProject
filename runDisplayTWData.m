% This code runs the main display program (displayTWData)

% The following programs should be downloaded and added to Matlab path 
% 1. Chronux: download from http://chronux.org/
% 2. ComonPrograms: https://github.com/supratimray/CommonPrograms
% 3. Add folder 'burstDurFiles' to path (Contains modified scripts of Gamma Length Project repository: https://github.com/supratimray/GammaLengthProjectCodes)
% 4. Supporting matlab files for MP: http://www.fuw.edu.pl/~durka/software/mp/.
% The data should be available in this folder

dataPath = 'G:\Bhargava_TravelingWaveProject\data\';

subjectName='alpaH'; expDate = '050817'; protocolName = 'GRF_002'; sizePos = 6; selectedElectrodes = [60 50 40 30 20 10];

%% Visualize time frequency and travelling wave plots 

displayTWData(subjectName,expDate,protocolName,dataPath,sizePos,selectedElectrodes);

%% Calculate Visualize travelling wave parameters
gridType = 'Microelectrode';
freqs = [30 60];
stimSize = sizePos;
[filtData,goodElectrodes,goodPos,timeVals] = loadLFPData(dataPath,subjectName,gridType,sizePos,freqs,1); %load lfp data
trialNum = 1; %calculating params for 1 trial
[outputs] = getTWCircParams(filtData(:,:,trialNum),timeVals,goodElectrodes,freqs,0); % for circ linear regression method
% [outputs] = getTWGradientParams(filtData(:,:,trialNum),goodElectrodes,timeVals,freqs,0); % for gradient method

%% Plot travelling wave parameters
sigPGD = outputs.pgd;
sigDirection = outputs.direction;
sigDirection(sigPGD<0.5) = nan; %consider only direction values corrosponding to significant pgd values (>0.5)
sigSpeed = outputs.speed;
sigSpeed(sigPGD<0.5) = nan; %consider only direction values corrosponding to significant pgd values (>0.5)
sigClusters = outputs.clusters;
sigClusters(sigPGD<0.5) = 0; %consider only electrode clusters corrosponding to significant pgd values (>0.5)

subplot(2,3,1:3)
plot(timeVals,sigPGD,'LineWidth',1,'Color','red')
hold on
yline(0.5)
xlabel('Time(s)')
yticks(0.5)
yticklabels('Significant PGD')
xlim([timeVals(1),timeVals(end)])
ylim([0 1])
hold off
title(['Phase Gradient Directionality for trial-',num2str(trialNum)])

subplot(2,3,4)
polarhistogram(sigDirection)
title('Direction distribution')

subplot(2,3,5)
histogram(sigSpeed,'Normalization','pdf')
xlabel('Speed (m/s)')
ylabel('PDF')
title('Speed distribution')

subplot(2,3,6)
plot(timeVals,sigClusters)
xlabel('Time(s)')
ylabel('Cluster Size')
title('TW ELectrode cluster at each time point') 