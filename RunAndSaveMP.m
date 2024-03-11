folderSourceString = 'E:\IISc_exp\IISC_work\GitScripts\monkeyData\'; % folder in which programs and data for this project are kept. Specify the string, e.g. 'E:\MonkeyData\TravellingWaveProject\';
%cd into the programs folder to run this script
%%%%%%%%%%%%%%%%%%%%%%%%%% Choice of parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Protocol details
monkeyName = 'alpaH'; expDate = '130418'; protocolName = 'GRF_004';
% monkeyName = 'kesariH'; expDate = '270218'; protocolName = 'GRF_001';
gridType = 'Microelectrode';
rfData = load([monkeyName 'MicroelectrodeRFData.mat']); % selecting good electrodes as per RMS values from Dubey and Ray, Sci Rep, 2020
goodElectrodes = rfData.highRMSElectrodes;

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

for i=1:numElectrodes
    eNum = goodElectrodes(i);

    if ~exist(fullfile(folderNameMP,(['elec' num2str(eNum) '/GaborMP/mp0.bok.000'])), 'file')

        tmp = load(fullfile(folderName,'segmentedData','LFP',['elec' num2str(eNum) '.mat']),'analogData');
        inputSignal = tmp.analogData;
    
        tag = sprintf('elec%d/',eNum);
    
        % Import the data
        X(:,:,1) = inputSignal';
        L = size(X,1);
        signalRange = [1 L]; % full range
        importData(X,folderNameMP,tag,signalRange,Fs);
    
        % perform Gabor decomposition
        Numb_points = L; % length of the signal
        runGabor(folderNameMP,tag,Numb_points,Max_iterations);

        gaborInfo = getGaborData(folderNameMP,tag,1);
        wrap=1;
        numTrials = length(gaborInfo);
        analogData = zeros(numTrials,L);
        
        for itrial=1:numTrials
            gaborData = gaborInfo{itrial}.gaborData;
            gaborData(:,(gaborData(2,:)==0)) = []; % Remove atoms with zero frequency
            analogData(itrial,:) = reconstructSignalFromAtomsMPP(gaborData,L,wrap);
        end

        filetoSave=fullfile(folderNameReconstMP,['elec' num2str(eNum) '.mat']);
        save(filetoSave,"analogData");
        clear gaborInfo
    else
        tag = sprintf('elec%d/',eNum);
        wrap = 1;
        gaborInfo = getGaborData(folderNameMP,tag,1);
        numTrials = length(gaborInfo);
        L = length(timeVals);        
        analogData = zeros(numTrials,L);
        
        for itrial=1:length(gaborInfo)
            gaborData = gaborInfo{itrial}.gaborData;
            gaborData(:,(gaborData(2,:)==0)) = []; % Remove atoms with zero frequency
            analogData(itrial,:) = reconstructSignalFromAtomsMPP(gaborData,L,wrap);
        end

        filetoSave=fullfile(folderNameReconstMP,['elec' num2str(eNum) '.mat']);
        save(filetoSave,"analogData");
        clear gaborInfo
    end
end
