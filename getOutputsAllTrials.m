dataPath = 'E:\IISc_exp\IISC_work\GitScripts\monkeyData\data\';
gridType = 'Microelectrode';
subjectName='alpaH'; expDate = '130418'; protocolName = 'GRF_004'; 
sizePos = 6;
req = 1; %load reconstructed MP data for whole gamma band

% get TW paramaters for whole gamma band
orientation = 1; % 0 (1)   22 (2)   45 (3)   67 (4)   90 (5)  112 (6)  135 (7)  157 (8) all orientations (9)
[allData,goodElectrodes,timeVals,~,~] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sizePos,orientation,1);
allData = permute(allData,[1,3,2]);

freqs = [25 75];
tic
for i = 1:size(allData,3)
    outputsG{i} = getTWCircParams(allData(:,:,i),timeVals,goodElectrodes,freqs,0,250);
end
toc
save('outputsWholeGamma','outputsG')
% get TW paramaters for slow gamma band
orientation = 1; % 0 (1)   22 (2)   45 (3)   67 (4)   90 (5)  112 (6)  135 (7)  157 (8) all orientations (9)
[allData,goodElectrodes,timeVals,rfData,parameters] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sizePos,orientation,2);
allData = permute(allData,[1,3,2]);

freqs = [25 45];
tic
for i = 1:size(allData,3)
    outputsSG{i} = getTWCircParams(allData(:,:,i),timeVals,goodElectrodes,freqs,0,250);
end
toc
save('outputsSlowGamma','outputsSG')


% get TW paramaters for fast gamma band
orientation = 1; % 0 (1)   22 (2)   45 (3)   67 (4)   90 (5)  112 (6)  135 (7)  157 (8) all orientations (9)
[allData,goodElectrodes,timeVals,rfData,parameters] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sizePos,orientation,3);
allData = permute(allData,[1,3,2]);

freqs = [45 75];
tic
for i = 1:size(allData,3)
    outputsHG{i} = getTWCircParams(allData(:,:,i),timeVals,goodElectrodes,freqs,0,250);
end
toc
save('outputsFastGamma','outputsHG')