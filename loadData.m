function [allData,goodElectrodes,timeVals,rfData,parameters] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sizePos,orientation,req)
    %Inputs
    %orientation = 0 deg (input 1)   22 deg (input 2)   45 deg (input 3)   67 deg (input 4)   90 deg (input 5)  112 deg (input 6)  135 deg (input 7)  157 deg (input 8) all orientations (input 9)
    %if req = 1, load the MP data of whole gamma band (25-75Hz)
    %req = 2, load the MP data of slow gamma band (25-45Hz)
    %req = 3, load the MP data of fast gamma band (45-75Hz)
    %req = 4, load unfiltered LFP data

    folderName = fullfile(dataPath,subjectName,gridType,expDate,protocolName);

    % Get good electrodes
    rfData = load([subjectName gridType 'RFData.mat']);
    goodElectrodes = rfData.highRMSElectrodes;
    goodElectrodes = goodElectrodes(goodElectrodes<=81); % Only microelectrodes
    numGoodElectrodes = length(goodElectrodes);

    % Get good trials
    parameters = load(fullfile(folderName,'extractedData','parameterCombinations.mat'));
    t = load(fullfile(folderName,'segmentedData','LFP','lfpInfo.mat'));
    timeVals = t.timeVals;

    badTrials = load(fullfile(folderName,'segmentedData','badTrials.mat'));
    badTrials = badTrials.badTrials;
    goodPos = setdiff(parameters.parameterCombinations{1,1,sizePos,1,orientation},badTrials);

    allData = zeros(numGoodElectrodes,length(goodPos),length(timeVals));

    if req==1
        dataReq = 'lfpMPGamma';
    elseif req==2
        dataReq = 'lfpMPSlowGamma';
    elseif req==3
        dataReq = 'lfpMPFastGamma';
    else
        dataReq = 'LFP';
    end
    for i=1:numGoodElectrodes
        lfpData = load(fullfile(folderName,'segmentedData',dataReq,['elec' num2str(goodElectrodes(i)) '.mat']));
        allData(i,:,:) = lfpData.analogData(goodPos,:);
    end
end