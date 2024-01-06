%overview of analysis
%load and filter data 
[filtData,goodElectrodes,goodPos,timeVals] = loadLFPData(dataPath,subjectName,gridType,stimSize,freqs);
%visualize the traveling wave
TWimagePalette(dataPath,subjectName,gridType,data,goodElectrodes,timeVals,timeFrame,videoTitle,trial,stimSize) 
%calculate TW parameters 
% set req to 1 for gradient method and to 2 for circ linear regression method
[outputs] = getTWParams(data,goodElectrodes,freqs,req,electrodeList);
%visualize the TW parameters
getTwparamsPlot(data,timeVals,goodElectrodes,times,outputs,sigValue,trial)