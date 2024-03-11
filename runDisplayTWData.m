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
sizePos = 7;
req = 1; %load reconstructed MP data for whole gamma band
selectedElectrodes = [40 41 42 43 44 45];
orientation = 9; % 0 (1)   22 (2)   45 (3)   67 (4)   90 (5)  112 (6)  135 (7)  157 (8) all orientations (9)
displayTWData(subjectName,expDate,protocolName,dataPath,sizePos,orientation,selectedElectrodes)

    
    