function plotStimLoc(dataPath,subjectName,gridType,sizeInd)
% using 150817/GRF_002 for this study presently
%     dataPath = 'G:\Bhargava_TravelingWaveProject\data\'; 
%     folderSourceString = 'G:\Bhargava_TravelingWaveProject\data\'; 
%     subjectName='alpaH'; gridType = 'Microelectrode';
    axisLims = [2 4.5 -4 -1.5];
    expDate = '050817'; 
    protocolName = 'GRF_002';
%     rfData = load('G:\Bhargava_TravelingWaveProject\data\AlpaH\Microelectrode\RFData.mat');
    rfData = load(fullfile(dataPath,subjectName,gridType,'\RFData.mat')); 
    
    %load visual loc for each electrode
    goodElectrodes = rfData.highRMSElectrodes;
    goodElectrodes = goodElectrodes(goodElectrodes<=81); % Only microelectrodes, ignoring ecog
    N = length(goodElectrodes);
    rfAziList = zeros(1,N);
    rfEleList = zeros(1,N);
        for i=1:N
            rfAziList(i) = rfData.rfStats(goodElectrodes(i)).meanAzi;
            rfEleList(i) = rfData.rfStats(goodElectrodes(i)).meanEle;
        end
        
    %load lfp data    
    folderName = fullfile(dataPath,subjectName,gridType,expDate,protocolName);
    parameters = load(fullfile(folderName,'extractedData','parameterCombinations.mat'));
    aValsUnique = parameters.aValsUnique;
    eValsUnique = parameters.eValsUnique;
    sValsUnique = 3*parameters.sValsUnique;
    plot(rfAziList,rfEleList,'+'); axis(axisLims); 
    hold on;
    plot(aValsUnique,eValsUnique,'ko');
    circle(aValsUnique,eValsUnique,sValsUnique(sizeInd));
    title('Stimulus location')
    xticks([2,2.5])
    xticklabels({'',''})
    yticks([-4,-3.5])
    yticklabels({'',''})
    axis square    
end
