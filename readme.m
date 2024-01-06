%overview of analysis
%% load and filter data 
dataPath = 'G:\Bhargava_TravelingWaveProject\data\';
subjectName='alpaH';
gridType = 'Microelectrode';
freqs =  [30 60];
req =1;
stimSize = 6;
[data,goodElectrodes,goodPos,timeVals] = loadLFPData(dataPath,subjectName,gridType,stimSize,freqs,req);

%% visualize the traveling wave
timeFrame = [0 1];
videoTitle = 'testTWVideo.avi';
trial = 1;
TWimagePalette(dataPath,subjectName,gridType,data,goodElectrodes,timeVals,timeFrame,videoTitle,trial,stimSize) 

%% calculate TW parameters 
% set req to 1 for gradient method and to 2 for circ linear regression method
req = 1;
electrodeList = [];
[outputs] = getTWParams(data,goodElectrodes,freqs,req,electrodeList);

%% visualize the TW parameters
sigValue = 0.5;
times = [0.055,0.06,0.065,0.07];
getTwparamsPlot(data,timeVals,goodElectrodes,times,outputs,sigValue,trial)
