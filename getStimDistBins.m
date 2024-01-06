function [goodDPos,distanceBins] = getStimDistBins(dataPath,subjectName,gridType,stimSize)
%load params
    expDate = '050817'; 
    protocolName = 'GRF_002';  
    folderName = fullfile(dataPath,subjectName,gridType,expDate,protocolName);
    rfData = load(fullfile(dataPath,subjectName,gridType,'\RFData.mat')); 
%% sort electrodes into distance bins
    goodElectrodes = rfData.highRMSElectrodes;
    goodElectrodes = goodElectrodes(goodElectrodes<=81); % Only microelectrodes
    N = length(goodElectrodes);
    rfAziList = zeros(1,N);
    rfEleList = zeros(1,N);
        for i=1:N
            rfAziList(i) = rfData.rfStats(goodElectrodes(i)).meanAzi;
            rfEleList(i) = rfData.rfStats(goodElectrodes(i)).meanEle;
        end
    % get the unique azi, ele and size info for the selected day
    parameters = load(fullfile(folderName,'extractedData','parameterCombinations.mat'));
    aValsUnique = parameters.aValsUnique;
    eValsUnique = parameters.eValsUnique;
    sValsUnique = 3*parameters.sValsUnique(stimSize);

    % Plot RFs and stimulus positions
    distanceList = sqrt((aValsUnique-rfAziList).^2+(eValsUnique-rfEleList).^2);

    distanceBins = 0:sValsUnique(1):1.5;
    NBins = length(distanceBins)-1;
    goodDPos = cell(1,NBins);
        for i=1:NBins
            goodDPos{i} = intersect(find(distanceList>distanceBins(i)),find(distanceList<=distanceBins(i+1)));
        end
end