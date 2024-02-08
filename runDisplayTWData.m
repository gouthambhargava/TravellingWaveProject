% This code runs the main display program (displayTWData)

% The following programs should be downloaded and added to Matlab path 
% 1. Chronux: download from http://chronux.org/
% 2. ComonPrograms: https://github.com/supratimray/CommonPrograms

% The data should be available in this folder
dataPath = 'D:\OneDrive - Indian Institute of Science\Supratim\Projects\Bhargava_TravelingWaveProject\data\';

subjectName='alpaH'; expDate = '050817'; protocolName = 'GRF_002'; sizePos = 6; selectedElectrodes = [60 50 40 30 20 10];

displayTWData(subjectName,expDate,protocolName,dataPath,sizePos,selectedElectrodes);
