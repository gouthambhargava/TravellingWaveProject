function [filtData,goodElectrodes,goodPos,timeVals] = loadLFPData(dataPath,subjectName,gridType,stimSize,freqs,req)
% load microelectrode data for specified  and filter lfp 
% using 150817/GRF_002 for this study presently
% inputs - 
%stimSize - give numerical value of size of stimulus being taken for the analysis (1-6)
% dataPath = 'G:\Bhargava_TravelingWaveProject\data\'; 
% subjectName='alpaH'; 
% gridType = 'Microelectrode'; 
% freqs - ex. [30 60];
% if req =1, get filtered data, else get raw data
    expDate = '050817'; 
    protocolName = 'GRF_002';  
    fs = 2000;
    folderName = fullfile(dataPath,subjectName,gridType,expDate,protocolName);
    rfData = load(fullfile(dataPath,subjectName,gridType,'\RFData.mat')); 
    goodElectrodes = rfData.highRMSElectrodes;
    goodElectrodes = goodElectrodes(goodElectrodes<=81); % Only microelectrodes
    
    % get the unique azi, ele and size info for the selected day
    parameters = load(fullfile(folderName,'extractedData','parameterCombinations.mat'));
    t = load(fullfile(folderName,'segmentedData','LFP','lfpInfo.mat'));
    timeVals = t.timeVals;
    badTrials = load(fullfile(folderName,'segmentedData','badTrials.mat'));
    badTrials = badTrials.badTrials;
    goodPos = parameters.parameterCombinations{1,1,stimSize,1,9};

    for i=1:length(goodElectrodes) 
        lfpData = load(fullfile(folderName,'segmentedData','LFP',['elec' num2str(goodElectrodes(i)) '.mat']));
        allTrials(:,:,i) = lfpData.analogData;
    end    
    allTrials = permute(allTrials,[3,2,1]);    

    data = allTrials(:,:,setdiff(goodPos,badTrials));
    if req~=1
        filtData = data;
    else
    normBand=freqs/(fs/2);
    filtOrder = 4;
    [b,a]=butter(filtOrder,normBand,'bandpass');    
    for i=1:size(data,3)
        signal = zscore(data(:,:,i))';
        filtData=filtfilt(b,a,signal)';
    end  
    end
end




