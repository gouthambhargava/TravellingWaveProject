% This code runs the main display program (displayTWData)
%% select data of particular orientation
dataPath = 'F:\monkeyData\data';
% dataPath = 'N:\Projects\Bhargava_TravelingWaveProject\data';
gridType = 'Microelectrode';

subjectName='alpaH'; expDate = '210817'; protocolName = 'GRF_002'; selectedElectrodes = [1 41 81]; freqRangeList{1} = [25 35]; freqRangeList{2} = [40 50];
stimPeriod = [0.25 0.75];
waveLengthLimit = 10;
sPos = 1; % spatial frequency: 0.5 (1), 1(2), 2 (3), 4 (4), 8 (5), all SFs (6). Note that the same code can be used for the size project also later where stimulus size is changed instead of spatial frequency
oriPos = 4; % orientation: 0 (1), 22.5 (2), 45 (3), 67.5 (4), 90 (5), 112.5 (6), 135 (7), 157.5 (8), all orientations (9)

analysisMethod = 'hilbert';
displayTWData(gridType,subjectName,expDate,protocolName,dataPath,sPos,oriPos,selectedElectrodes,analysisMethod,freqRangeList,waveLengthLimit,loadDataFlag);