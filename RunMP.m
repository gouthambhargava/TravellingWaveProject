folderSourceString = 'E:\IISc_exp\IISC_work\GitScripts\monkeyData\'; % folder in which programs and data for this project are kept. Specify the string, e.g. 'E:\MonkeyData\TravellingWaveProject\';
%cd into the programs folder to run this script
%%%%%%%%%%%%%%%%%%%%%%%%%% Choice of parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Protocol details
monkeyName = 'alpaH'; expDate = '130418'; protocolName = 'GRF_004';
% monkeyName = 'kesariH'; expDate = '270218'; protocolName = 'GRF_001';
gridType = 'Microelectrode';
rfData = load([monkeyName 'MicroelectrodeRFData.mat']); % selecting good electrodes as per RMS values from Dubey and Ray, Sci Rep, 2020
goodElectrodes = rfData.highRMSElectrodes;
freqRange = [30 60];

numElectrodes = length(goodElectrodes);

folderName = fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName);
folderNameMP = fullfile(folderSourceString,'/data',monkeyName,gridType,expDate,protocolName,'mpAnalysis'); % The path must be relative when using Windows (i.e., should not start with c: etc. So this part only works when we are in the programs folder)
makeDirectory(folderNameMP);
folderNameReconstMP = fullfile(folderSourceString,'/data',monkeyName,gridType,expDate,protocolName,'segmentedData','lfpMP');
makeDirectory(folderNameReconstMP);

tmp = load(fullfile(folderName,'segmentedData','LFP','lfpInfo.mat'),'timeVals');
timeVals = tmp.timeVals;
Fs = round(1/((timeVals(2)-timeVals(1))));
Max_iterations = 500; % number of iterations
adaptiveDictionaryParam=0.9;
dictionarySize=2500000;

%% run this section to generate gaborData and save it
for i=1:numElectrodes
    eNum = goodElectrodes(i);

    tmp = load(fullfile(folderName,'segmentedData','LFP',['elec' num2str(eNum) '.mat']),'analogData');
    inputSignal = tmp.analogData;
    tag = sprintf('elec%d/',eNum);
    
    % perform Gabor decomposition
    [gaborInfo,header] = getStochasticDictionaryMP3p1(inputSignal,timeVals,maxIteration,adaptiveDictionaryParam,dictionarySize);
    filetoSave=fullfile(folderNameMP,['elec' num2str(eNum) '.mat']);
    save(filetoSave,'gaborInfo','header');
    clear gaborInfo
end 

%% run this section if MP data has been saved and LFP time series is to be obtained from the MP data
% modify freqRange accordingly to obtain time series comprising of only
% those atoms
for i=1:numElectrodes
    eNum = goodElectrodes(i);
    fileName=fullfile(folderNameReconstMP,['elec' num2str(eNum) '.mat']);
    load(filetoSave)
    % some constants related to the internal structure of the 'book':
    scale=1;freq=2;pos=3;mod=4;amp=5;phi=6;HSr=1;HSs=2;
    reconSig=zeros(1,h(HSs));
    numTrials = size(gaborInfo,1);
    
    for ind=1:numTrials
        gaborInt = squeeze(gaborInfo(ind,:,:));
        freqValues = cat(1,find(gaborInt(:,2)<freqRange(1)),find(gaborInt(:,2)>freqRange(2)));
        gaborInt(freqValues,:) = [];
        reconSig = reconSig+gabor(header(1,HSs)/header(1,HSr),header(1,HSr),gaborInt(:,scale),gaborInt(:,freq),gaborInt(:,pos),gaborInt(:,amp),gaborInt(:,phi));
    end
    filetoSave=fullfile(folderNameReconstMP,['LFP_elec' num2str(eNum) '.mat']);
    save(filetoSave,'reconSig');
end
    