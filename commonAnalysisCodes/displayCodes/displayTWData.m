% Analysis Method

% Four analyses are done: 1) TF spectrum, 2) filtering the signals in
% prespecified frequency bands, 3) finding epochs where bursts are detected
% and 4) traveling wave characteristics

% "hilbert" method involves using Multi-taper, traditional filtering and Hilbert transform for analyses 1-3
% "Matching Pursuit" method involves using MP for TF spectrum, reconstructing the
% signal within a pre-specified band by using atoms with centre frequencies
% within tha band, and doing burst estimation using MP as well.

% Burst characteristics need to be calculated for all trials
% simultaneously. To find common time epochs where traveling waves are
% found, this analysis needs to be done for all electrodes also. These
% analyses are therefore done upfront even though results are shown only
% for some electrodes and trials.

function displayTWData(gridType,subjectName,expDate,protocolName,dataPath,sPos,oriPos,selectedElectrodes,analysisMethod,freqRangeList,lengthLimit)

fontSizeSmall = 10; fontSizeMedium = 12; fontSizeLarge = 16;fontSizeTiny = 6;
backgroundColor = 'w'; panelHeight = 0.125;


numFrequencyRanges = length(freqRangeList);
colorNamesFreqRanges = gray(numFrequencyRanges+1);
% waveColors = cat(1,[0 0 1],[0 1 0],[0 1 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(gridType,'Microelectrode')
    disp('Getting data...');
    [allData,goodElectrodes,timeVals,rfData,parameters] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sPos,oriPos);
    % else
    %     % load EEG data
end
numGoodElectrodes = length(goodElectrodes);

%%%%%%%%%%%%%%%%%%%%%%% Do burst estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change - burst estimation has been moved to a different script to
% accomodate other methods of bursts detection as well. Bursts will be
% detected for single trials only, so they have been moved to
% plot_callback1
%%%%%%%%%%%%%%%%%%%%%%%% Plot RF information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hGridPlots = getPlotHandles(1,2,[0.025 0.65 0.3 0.2],0.01,0.02,0);
numSelectedElectrodes = length(selectedElectrodes);
colorNamesElectrodes = jet(numSelectedElectrodes);
stimulusDurationS = [0 0.8]; % Stimulus duration to be highlighted
baselinePeriodS = [-0.5 0];
stimulusPeriodS = [0.25 0.75];
analysisPeriodS = [-0.5 1];
thresholdFactor = 3;
filterOrder = 4;
electrodeArray = showRFPositionsSelectedElectrodes(hGridPlots,goodElectrodes,selectedElectrodes,rfData,parameters,colorNamesElectrodes);

% get location for each electrode for segmentation into clusters
locList = nan(numGoodElectrodes,2);
for iElec = 1:numGoodElectrodes
    [locList(iElec,1),locList(iElec,2)] = find(electrodeArray==goodElectrodes(iElec));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%  Data selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
hPanel1 = uipanel('Title','Data','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0 1-panelHeight 0.18 panelHeight]);

% Selected electrodes
uicontrol('Parent',hPanel1,'Unit','Normalized','Position',[0 0.5 0.25 0.5],'Style','text','String','Elecs','FontSize',fontSizeMedium);
hElectrodes = uicontrol('Parent',hPanel1,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.25 0.5 0.75 0.5],'Style','edit','String',num2str(selectedElectrodes),'FontSize',fontSizeMedium);

% Trial Number
numTrials = size(allData,2);
uicontrol('Parent',hPanel1,'Unit','Normalized','Position',[0 0 0.25 0.5],'Style','text','String','Trial','FontSize',fontSizeMedium);
trialNumList = 1:numTrials;
hTrial = uicontrol('Parent',hPanel1,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.25 0 0.15 0.5],'Style','popup','String',trialNumList,'FontSize',fontSizeMedium);

% Plot data
uicontrol('Parent',hPanel1,'Unit','Normalized','Position',[0.4 0 0.2 0.5],'Style','pushbutton','String','plot','FontSize',fontSizeSmall,'Callback',{@plot_Callback});
uicontrol('Parent',hPanel1,'Unit','Normalized','Position',[0.6 0 0.2 0.5],'Style','pushbutton','String','Rescale','FontSize',fontSizeSmall,'Callback',{@rescale_Callback});
uicontrol('Parent',hPanel1,'Unit','Normalized','Position',[0.8 0 0.2 0.5],'Style','pushbutton','String','Clear','FontSize',fontSizeSmall,'Callback',{@cla_Callback});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Axis Ranges %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPanel2 = uipanel('Title','AxisRanges1','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.18 1-panelHeight 0.15 panelHeight]);
axisRange1List0{1} = [0 80]; axisRange1Name{1} = 'FreqLims (Hz)';
axisRange1List0{2} = [-0.5 1]; axisRange1Name{2} = 'TimeLims (S)';
axisRange1List0{3} = [-2 2]; axisRange1Name{3} = 'cLims (raw)';

numAxisRanges1 = length(axisRange1List0);
hAxisRange1Min = cell(1,numAxisRanges1);
hAxisRange1Max = cell(1,numAxisRanges1);

for ii=1:numAxisRanges1
    uicontrol('Parent',hPanel2,'Unit','Normalized','Position',[0 1-ii/numAxisRanges1 0.5 1/numAxisRanges1],'Style','text','String',axisRange1Name{ii},'FontSize',fontSizeSmall);
    hAxisRange1Min{ii} = uicontrol('Parent',hPanel2,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.5 1-ii/numAxisRanges1 0.25 1/numAxisRanges1], ...
        'Style','edit','String',num2str(axisRange1List0{ii}(1)),'FontSize',fontSizeSmall);
    hAxisRange1Max{ii} = uicontrol('Parent',hPanel2,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.75 1-ii/numAxisRanges1 0.25 1/numAxisRanges1], ...
        'Style','edit','String',num2str(axisRange1List0{ii}(2)),'FontSize',fontSizeSmall);
end

hPanel3 = uipanel('Title','AxisRanges2','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.33 1-panelHeight 0.16 panelHeight]);

axisRange2List0 = cell(1,numFrequencyRanges);
axisRange2Name = cell(1,numFrequencyRanges);
for iFreq=1:numFrequencyRanges
    axisRange2List0{iFreq} = [-50 50];
    axisRange2Name{iFreq} = ['yRange (' num2str(freqRangeList{iFreq}(1)) '-' num2str(freqRangeList{iFreq}(2)) ' Hz)'];
end

hAxisRange2Min = cell(1,numFrequencyRanges);
hAxisRange2Max = cell(1,numFrequencyRanges);

for ii=1:numFrequencyRanges
    uicontrol('Parent',hPanel3,'Unit','Normalized','Position',[0 1-ii/numFrequencyRanges 0.5 1/numFrequencyRanges],'Style','text','String',axisRange2Name{ii},'FontSize',fontSizeSmall);
    hAxisRange2Min{ii} = uicontrol('Parent',hPanel3,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.5 1-ii/numFrequencyRanges 0.25 1/numFrequencyRanges], ...
        'Style','edit','String',num2str(axisRange2List0{ii}(1)),'FontSize',fontSizeSmall);
    hAxisRange2Max{ii} = uicontrol('Parent',hPanel3,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.75 1-ii/numFrequencyRanges 0.25 1/numFrequencyRanges], ...
        'Style','edit','String',num2str(axisRange2List0{ii}(2)),'FontSize',fontSizeSmall);
end

%%%%%%%%%%%%%%%%%%%%  Traveling Wave plot panel %%%%%%%%%%%%%%%%%%%%%%%%

hPanel4 = uipanel('Title','Traveling Wave','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.49 1-panelHeight 0.31 panelHeight]);

uicontrol('Parent',hPanel4,'Unit','Normalized','Position',[0 0.5 0.2 0.5],'Style','text','String','Electrode Fraction','FontSize',fontSizeSmall);
hElecFrac = uicontrol('Parent',hPanel4,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.2 0.7 0.15 0.3], ...
    'Style','edit','String','0.5','FontSize',fontSizeSmall);

electrodeChoices = {'selected','all'};
uicontrol('Parent',hPanel4,'Unit','Normalized','Position',[0 0 0.2 0.5],'Style','text','String','Select Electrodes','FontSize',fontSizeSmall);
hSelectElec = uicontrol('Parent',hPanel4,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.2 0.2 0.15 0.3], ...
    'Style','popup','String',electrodeChoices,'FontSize',fontSizeSmall);

waveSegOptions = {'Simple Segmentation','Point-point wobble','Full segment wobble','Wave strength'};
% note - wave strength is to be implemented for cluster wise PGD/direction
% calculation - for EEG
uicontrol('Parent',hPanel4,'Unit','Normalized','Position',[0.37 0.55 0.15 0.5],'Style','text','String','Wave Seg Option','FontSize',fontSizeSmall);
hSegReq = uicontrol('Parent',hPanel4,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.52 0.7 0.15 0.3], ...
    'Style','popup','String',waveSegOptions,'FontSize',fontSizeSmall);

waveMethods = {'Simple circLin regression','Cluster circLin regression'};
% current version only performs simple circular linear regression - cluster
% based circlinear regression to be implemented for EEG.
uicontrol('Parent',hPanel4,'Unit','Normalized','Position',[0.37 0 0.15 0.5],'Style','text','String','Wave Detection','FontSize',fontSizeSmall);
hWaveMethod = uicontrol('Parent',hPanel4,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.52 0.2 0.15 0.3], ...
    'Style','popup','String',waveMethods,'FontSize',fontSizeSmall);

uicontrol('Parent',hPanel4,'Unit','Normalized','Position',[0.68 0.55 0.15 0.5],'Style','text','String','Wave wobble','FontSize',fontSizeSmall);
hWobbleReq = uicontrol('Parent',hPanel4,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.83 0.7 0.15 0.3], ...
    'Style','edit','String','5','FontSize',fontSizeSmall);

uicontrol('Parent',hPanel4,'Unit','Normalized','Position',[0.68 0 0.15 0.5],'Style','text','String','Ref Choice','FontSize',fontSizeSmall);
hSelectRef = uicontrol('Parent',hPanel4,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.83 0.2 0.15 0.3], ...
    'Style','popup','String',{'Avg',num2str(goodElectrodes')},'FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%  Phase propagation plot panel %%%%%%%%%%%%%%%%%%%%%%%%
hPanel5 = uipanel('Title','Plot Phase Propagation','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.8 1-panelHeight 0.2 panelHeight]);

timeRangeprop0 = [0 1];
uicontrol('Parent',hPanel5,'Unit','Normalized','Position',[0 0.5 0.5 0.5],'Style','text','String','Time Range (s)','FontSize',fontSizeMedium);
hTimeRangePropMin = uicontrol('Parent',hPanel5,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.5 0.5 0.25 0.5], ...
    'Style','edit','String',num2str(timeRangeprop0(1)),'FontSize',fontSizeSmall);
hTimeRangePropMax = uicontrol('Parent',hPanel5,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.75 0.5 0.25 0.5], ...
    'Style','edit','String',num2str(timeRangeprop0(2)),'FontSize',fontSizeSmall);

plotToggle = uicontrol('Parent',hPanel5,'Unit','Normalized','Position',[0 0 0.5 0.5],'Style','togglebutton','String','Plot/Pause','FontSize',fontSizeMedium,'Callback',{@plot_Callback2},'Value',0);
uicontrol('Parent',hPanel5,'Unit','Normalized','Position',[0.5 0 0.5 0.5],'Style','pushbutton','String','Clear','FontSize',fontSizeMedium,'Callback',{@cla_Callback2});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot handles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hTF = getPlotHandles(numSelectedElectrodes+1,2,[0.025 0.05 0.3 0.55],0.025,0.025);
hStats = getPlotHandles(4,numFrequencyRanges,[0.35 0.52 0.425 0.33],0.025,0.01);
hBursts = getPlotHandles(1,numFrequencyRanges,[0.35 0.33 0.425 0.17],0.025,0.01);
hSignal = getPlotHandles(numSelectedElectrodes,numFrequencyRanges,[0.35 0.05 0.425 0.25],0.025,0.01);
hGridPlots2 = getPlotHandles(3,numFrequencyRanges,[0.8 0.05 0.19 0.8]);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Needs to be plotted only once
plotTFMT(hTF(1,1),allData,timeVals,axisRange1List0,freqRangeList,colorNamesFreqRanges,'raw',[],stimulusDurationS); % Time-frequency power spectrum for all trials

outputsTW = cell(1,numFrequencyRanges);
state = 0; % state variable to display the phases
waveVector = cell(1,numFrequencyRanges);
uniqueDirs = cell(1,numFrequencyRanges);
% lengthLimit = 10;


    function plot_Callback(~,~)

        % Data choices
        selectedElectrodes0 = str2num(get(hElectrodes,'String')); %#ok<ST2NM>
        trialNum = trialNumList((get(hTrial,'val')));
        wobbleLim = str2double(get(hWobbleReq,'String'));
        segOption = get(hSegReq,'val');
        waveMethod = get(hWaveMethod,'val');
        refPhaseChoice = get(hSelectRef,'val');
        gridLayout = rot90(reshape(1:81,[9,9]),2); %set the grid layout, change to add EEG 

        % define color scheme for waves
        waveColors = parula(3); % 1- planar 2 -spiral 3 - complex
        allWaveColors = cell(1,numFrequencyRanges);

        disp('Performing burst analysis...');
        filteredSignal = zeros(numGoodElectrodes,size(allData,2),size(allData,3),numFrequencyRanges);
        burstTS = nan(numGoodElectrodes,size(allData,2),size(allData,3),numFrequencyRanges);
        phaseMatrix = nan(numGoodElectrodes,size(allData,3),numFrequencyRanges);
        refCube = nan(size(gridLayout,1),size(gridLayout,2),size(allData,3),numFrequencyRanges);
        burstMatrix = nan(numGoodElectrodes,size(allData,3),numFrequencyRanges);
        
        for i=1:numFrequencyRanges
            for j=1:numGoodElectrodes
                [~,~,~,burstTS(j,:,:,i),filteredSignal(j,:,:,i),~] = getHilbertBurst(squeeze(allData(j,:,:)),timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,freqRangeList{i},filterOrder,1,analysisPeriodS);
            end
            burstMatrix(:,:,i) = burstTS(:,trialNum,:,i);
            phaseMatrix(:,:,i) = angle(hilbert(squeeze(filteredSignal(:,trialNum,:,i))'))';
            if refPhaseChoice==1
                refPhase = circMeanNan(phaseMatrix(:,:,i),1);
            else
                refPhase = phaseMatrix(redPhaseChoice-1,:,i);
            end
             refMatrix = angle(exp(1i*(phaseMatrix(:,:,i) - refPhase)));
            for e=1:numGoodElectrodes
                refCube(locList(e,1),locList(e,2),:,i) = refMatrix(e,:);
            end    
        end
        filteredSignalTF = squeeze(filteredSignal(:,trialNum,:,:));

        %get burst fraction statistics
        elecFrac = str2double(get(hElecFrac,'string'));
        elecChoice = electrodeChoices{get(hSelectElec,'val')};

        %% Calculate TW parameters
        % the GUI now takes analysed data to show plots
        % as waveMethod 2 takes 10-15 minutes per trial. This has been
        % applied for wave method 1 as well.
        % naming convention for data - subject name/SFPos/oriPos/electrode fraction/electrode selection/wave method
        % ex alpaH_24_0.5T_elecChoice1_met2
        disp('Getting TW parameters');
        % waveVector = nan(3,length(timeVals),numFrequencyRanges);
        filename = [subjectName,'_',num2str(oriPos),num2str(sPos),'_',num2str(elecFrac),'T_',elecChoice,'_met',num2str(waveMethod),'.mat'];
        if exist(filename,'file')>0
        outputs = load(filename);
            outputsTW = outputs.outputs(:,trialNum);
            for i = 1:numFrequencyRanges
                burstVec = nansum(burstMatrix(:,:,i))/numGoodElectrodes;
                burstVec(burstVec>=elecFrac) = 1;
                burstVec(burstVec~=1) = nan;
                outputsTW{i}.burstVec = burstVec;
            end
        else
            for j = 1:numTrials
                for i = 1:numFrequencyRanges
                    bursts = squeeze(burstTS(:,j,:,i));
                    phases = angle(hilbert(squeeze(filteredSignal(:,j,:,i))'))';
                outputs{i,j} = getTWCircParams(phases,bursts,timeVals,locList,elecFrac,elecChoice,gridType,waveMethod);
                end
            end
            cd('dualGammaWaveProject/data')
            save(filename,'outputs')
            outputsTW = outputs{:,trialNum};
        end

        for i = 1:numFrequencyRanges
            [waveVector{i},uniqueDirs{i},~,clusterCube,waveType,directionCube,phaseCube] = getWaveParameters(outputsTW{i},phaseMatrix(:,:,i),locList,timeVals,lengthLimit, wobbleLim, segOption,stimulusPeriodS,waveMethod);
            if isempty(phaseCube)
                phaseCube = nan(size(gridLayout,1),size(gridLayout,1),length(timeVals));
                clusterCube = nan(size(gridLayout,1),size(gridLayout,1),length(timeVals));
                directionCube = nan(size(gridLayout,1),size(gridLayout,1),length(timeVals));
            end   
            outputsTW{i}.clusterCube = clusterCube;
            outputsTW{i}.directionCube = directionCube;
            outputsTW{i}.phaseCube = phaseCube;  
            outputsTW{i}.refCube = refCube(:,:,:,i);
            waveColorsInt = nan(numel(waveType),3);            
            for j = 1:numel(waveType)
                waveColorsInt(j,:) = waveColors(waveType(j),:);
            end
            allWaveColors{i} = waveColorsInt;
        end
        disp('TW analysis completed')
        %%
        
        if ~isequal(selectedElectrodes0,selectedElectrodes)
            showRFPositionsSelectedElectrodes(hGridPlots(1,:),goodElectrodes,selectedElectrodes0,rfData,parameters,colorNamesElectrodes);
            selectedElectrodes = selectedElectrodes0;
        end

        % Axis Ranges
        axisRange1List = cell(1,numAxisRanges1);
        for i=1:numAxisRanges1
            axisRange1List{i} = [str2double(get(hAxisRange1Min{i},'String')) str2double(get(hAxisRange1Max{i},'String'))];
        end
        axisRange2List = cell(1,numFrequencyRanges);
        for i=1:numFrequencyRanges
            axisRange2List{i} = [str2double(get(hAxisRange2Min{i},'String')) str2double(get(hAxisRange2Max{i},'String'))];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%% Plot data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plot TF data for all channels
        plotTFMT(hTF(1,2),allData,timeVals,axisRange1List,freqRangeList,colorNamesFreqRanges,'raw',trialNum,stimulusDurationS);

        %Plot TF info for selected channels
        for i=1:numSelectedElectrodes
            ePos = find(selectedElectrodes(i)==goodElectrodes);
            signalAllTrials = allData(ePos,:,:);

            if strcmp(analysisMethod,'hilbert') % Use Multi-taper method
                plotTFMT(hTF(i+1,1),signalAllTrials,timeVals,axisRange1List,freqRangeList,colorNamesFreqRanges,'raw',[],stimulusDurationS); % Time-frequency power spectrum for all trials
                plotTFMT(hTF(i+1,2),signalAllTrials,timeVals,axisRange1List,freqRangeList,colorNamesFreqRanges,'raw',trialNum,stimulusDurationS);
            else
                plotTFMP(hTF(i+1,1),subjectName,expDate,protocolName,dataPath,gridType,sPos,oriPos,selectedElectrodes(i),timeVals,axisRange1List,freqRangeList,colorNamesFreqRanges); % Write this code - this data should be saved somewhere
            end
            hold(hTF(i+1,2),'on');

            % Plot filtered signal and show bursts
            for j=1:numFrequencyRanges
                plot(hSignal(i,j),timeVals,squeeze(filteredSignalTF(ePos,:,j)),'color',colorNamesFreqRanges(j,:));
                hold(hSignal(i,j),'on');
                plot(hSignal(i,j),timeVals,squeeze(burstMatrix(ePos,:,j))-1,'color','r','linewidth',2);
                axis(hSignal(i,j),[axisRange1List{2} axisRange2List{j}]);
                plot(hTF(i+1,2),timeVals, mean(freqRangeList{j})+ burstMatrix(ePos,:,j)-1,'color',colorNamesFreqRanges(j,:),'linewidth',2);
            end
        end

        % Plot TW stats
        for i = 1:numFrequencyRanges
            % Plot bursts
            chanVals = repmat((1:numel(goodElectrodes))',1,size(burstMatrix,2));
            plot(timeVals,burstMatrix(:,:,i)+chanVals-1,'parent',hBursts(1,i),'color','r','linewidth',1);
            hold(hBursts(1,i),'on');

            % Show time points where the grid as a whole has bursts
            plot(hBursts(1,i),timeVals(outputsTW{i}.burstVec==1),numel(goodElectrodes)+1,'|','Color','black');
            xlim(hBursts(1,i),axisRange1List{2});

            % plot different wave segments with the bursts - colored individually
            for k = 1:size(waveVector{i},1)
                plot(hBursts(1,i),timeVals,squeeze(waveVector{i}(k,:))*0,'|','Color',allWaveColors{i}(k,:),'linewidth',1);
                hold(hBursts(1,i),'on');
            end

            % plot the rest of the TW parameters
            pgdTemp = outputsTW{i}.pgd;
            if min(size(pgdTemp)) == 1
                pgdTemp = reshape(pgdTemp,[1,length(pgdTemp)]);
            else
                pgdTemp = nanmean(pgdTemp);
            end
            pgdTemp(isnan(pgdTemp)) = 0;
            plot(hStats(1,i),timeVals,pgdTemp,'color',colorNamesFreqRanges(i,:));
            axis(hStats(1,i),[axisRange1List{2},0 1]);
            xticks(hStats(1,i),0);
            xticklabels(hStats(1,i),'');
            legend(hStats(1,i),{'PGD'},'Location','best','AutoUpdate','off')

            speedTemp = outputsTW{i}.speed;
            speedTemp(isnan(speedTemp)) = 0;
            if size(speedTemp,1)>1
                speedTemp = mean(speedTemp);
            end
            plot(hStats(2,i),timeVals,speedTemp,'color',colorNamesFreqRanges(i,:));
            axis(hStats(2,i),[axisRange1List{2},0 1]);
            xticks(hStats(2,i),0);
            xticklabels(hStats(2,i),'');
            legend(hStats(2,i),{'Speed'},'Location','best','AutoUpdate','off')

            dirTemp = wrapToPi(outputsTW{i}.direction);
            if min(size(dirTemp))>1
                dirTemp = circMeanNan(dirTemp);
            end
            plot(hStats(3,i),timeVals,dirTemp,'color',colorNamesFreqRanges(i,:));
            axis(hStats(3,i),[axisRange1List{2},-4 4]);
            xticks(hStats(3,i),0);
            xticklabels(hStats(3,i),'');
            legend(hStats(3,i),{'Direction'},'Location','best','AutoUpdate','off')

            plot(hStats(4,i),timeVals,outputsTW{i}.coh,'color','r');
            hold(hStats(4,i),'on');
            plot(hStats(4,i),timeVals,sum(burstMatrix(:,:,i),'omitnan')/numGoodElectrodes,'color','k');
            line(axisRange1List{2},[elecFrac elecFrac],'parent',hStats(4,i),'color','b');
            axis(hStats(4,i),[axisRange1List{2} 0 1]);
            legend(hStats(4,i),{'Coh','Bursts'},'Location','best','AutoUpdate','off')

            % Plot angle plots for calculated directions
            makePolarPlot(uniqueDirs{i},10,hGridPlots2(3,i),allWaveColors{i},1)
        end

        htextPanel1 = uipanel('Unit','Normalized','Position',[0.89 0.05 0.02 0.02]);
        uicontrol('Parent',htextPanel1,'Unit','Normalized','Position',[0 0 1 1],'Style','text','String',num2str(['270',char(176)]),'FontSize',fontSizeTiny);

        htextPanel2 = uipanel('Unit','Normalized','Position',[0.89 0.27 0.02 0.02]);
        uicontrol('Parent',htextPanel2,'Unit','Normalized','Position',[0 0 1 1],'Style','text','String',num2str(['90',char(176)]),'FontSize',fontSizeTiny);

        htextPanel3 = uipanel('Unit','Normalized','Position',[0.78 0.16 0.02 0.02]);
        uicontrol('Parent',htextPanel3,'Unit','Normalized','Position',[0 0 1 1],'Style','text','String',num2str(['180',char(176)]),'FontSize',fontSizeTiny);

        htextPanel4 = uipanel('Unit','Normalized','Position',[0.99 0.16 0.02 0.02]);
        uicontrol('Parent',htextPanel4,'Unit','Normalized','Position',[0 0 1 1],'Style','text','String',num2str(['0',char(176)]),'FontSize',fontSizeTiny);


    end
%%
    function plot_Callback2(~,~)

        % get inputs
        waveMethod = get(hWaveMethod,'val'); % to be implemented in future versions
        [X,Y] = meshgrid(1:1:9);

        if get(plotToggle,'Value') == 1
            state = 1; % This allows the program to start displaying the phase plots
            timeRangeProp = [str2double(get(hTimeRangePropMin,'String')) str2double(get(hTimeRangePropMax,'String'))];
            timeRangeProp = dsearchn(timeVals',timeRangeProp(1)):dsearchn(timeVals',timeRangeProp(2));

            %%%%%%%%%%%% Plot single trial and indicate bursts %%%%%%%%%%%%
            % Initialize
            for j=1:numFrequencyRanges
                xLinesBurst(1,j) = line([timeVals(timeRangeProp(1)) timeVals(timeRangeProp(1))],get(hBursts(1,j),'YLim'),'parent',hBursts(1,j),'LineWidth',2,'Color','black'); %#ok<*AGROW>
                for k=1:size(hSignal,1)
                    xLinesSignal(k,j) = line([timeVals(timeRangeProp(1)) timeVals(timeRangeProp(1))],get(hSignal(k,j),'YLim'),'parent',hSignal(k,j),'LineWidth',2,'Color','black'); %#ok<*AGROW>
                end
                for k=1:size(hStats,1)
                    xLinesStats(k,j) = line([timeVals(timeRangeProp(1)) timeVals(timeRangeProp(1))],get(hStats(k,j),'YLim'),'parent',hStats(k,j),'LineWidth',2,'Color','black');
                end
            end

            for ind = 1:numel(timeRangeProp)
                if state==0 % if state changes to 0, delete lines and return
                    for j=1:numFrequencyRanges
                        delete(xLinesBurst(1,j));
                        for k=1:size(hSignal,1)
                            delete(xLinesSignal(k,j));
                        end
                        for k=1:size(hStats,1)
                            delete(xLinesStats(k,j));
                        end
                    end
                    return;
                end

                for j=1:numFrequencyRanges
                    directionCube = outputsTW{j}.directionCube;
                    phaseCube = outputsTW{j}.phaseCube;
                    if waveMethod==1 % uncomment to plot the relative
                    % phase for method 1
                        propPlotSec = cos(outputsTW{j}.refCube);
                        plotLims = [-1 1];
                    else
                        propPlotSec = outputsTW{j}.clusterCube;
                        plotLims = [0 1];
                    end
                    
                    % Delete old lines and create new ones
                    delete(xLinesBurst(1,j));
                    xLinesBurst(1,j) = line([timeVals(timeRangeProp(ind)) timeVals(timeRangeProp(ind))],get(hBursts(1,j),'YLim'),'parent',hBursts(1,j),'LineWidth',2,'Color','black');

                    for k=1:size(hSignal,1)
                        delete(xLinesSignal(k,j));
                        xLinesSignal(k,j) = line([timeVals(timeRangeProp(ind)) timeVals(timeRangeProp(ind))],get(hSignal(k,j),'YLim'),'parent',hSignal(k,j),'LineWidth',2,'Color','black');
                    end
                    for k=1:size(hStats,1)
                        delete(xLinesStats(k,j));
                        xLinesStats(k,j) = line([timeVals(timeRangeProp(ind)) timeVals(timeRangeProp(ind))],get(hStats(k,j),'YLim'),'parent',hStats(k,j),'LineWidth',2,'Color','black');
                    end

                    % plot absolute phases
                    imagesc(cos(phaseCube(:,:,timeRangeProp(ind))),'parent',hGridPlots2(1,j));
                    colormap(hGridPlots2(1,j),'parula');
                    clim(hGridPlots2(1,j),[-1 1]);
                    axis(hGridPlots2(1,j),'off')
                    % colorbar(hGridPlots2(1,j),'northoutside');
                    hold(hGridPlots2(1,j),'on')
                    directionGrid = directionCube(:,:,timeRangeProp(ind));
                    quiver(X,Y,cos(directionGrid),sin(directionGrid),'Color','white','LineWidth',1,'AutoScaleFactor',0.5,'parent',hGridPlots2(1,j))
                    axis(hGridPlots2(1,j),'off')
                    axis(hGridPlots2(1,j),'square')
                    set(hGridPlots2(1,j),'YDir','Normal')

                    if waveMethod==1 %uncomment to plot the relative     
                        imagesc(propPlotSec(:,:,timeRangeProp(ind)),'parent',hGridPlots2(2,j));
                        set(hGridPlots2(2,j),'YDir','Normal')
                        quiver(X,Y,cos(propPlotSec(:,:,timeRangeProp(ind))),sin(propPlotSec(:,:,timeRangeProp(ind))),'Color','white','LineWidth',1,'AutoScaleFactor',0.5,'parent',hGridPlots2(2,j))
                    else
                        pcolor(propPlotSec(:,:,timeRangeProp(ind)),'parent',hGridPlots2(2,j));
                    end
                    colormap(hGridPlots2(2,j),'jet');
                    clim(hGridPlots2(2,j),plotLims);
                    hold(hGridPlots2(2,j),'on')
                    axis(hGridPlots2(2,j),'off')
                    axis(hGridPlots2(2,j),'square')
                    colorbar(hGridPlots2(2,j),'northoutside');
                end
                drawnow
            end
        else
            state = 0;
        end
    end

    function cla_Callback(~,~)
        claGivenPlotHandle(hTF(2:end,:)); cla(hTF(1,2));
        claGivenPlotHandle(hSignal);
        claGivenPlotHandle(hStats);
        claGivenPlotHandle(hGridPlots2([3,6]))
        claGivenPlotHandle(hBursts)

    end
    function rescale_Callback(~,~)
        freqLims = [str2double(get(hAxisRange1Min{1},'String')) str2double(get(hAxisRange1Max{1},'String'))];
        timeLims = [str2double(get(hAxisRange1Min{2},'String')) str2double(get(hAxisRange1Max{2},'String'))];
        cLims = [str2double(get(hAxisRange1Min{3},'String')) str2double(get(hAxisRange1Max{3},'String'))];

        % Rescale TF plots
        rescaleGivenPlotHandle(hTF,[timeLims freqLims]);
        rescaleZGivenPlotHandle(hTF,cLims);

        for i=1:numFrequencyRanges
            yLims = [str2double(get(hAxisRange2Min{i},'String')) str2double(get(hAxisRange2Max{i},'String'))];
            rescaleGivenPlotHandle(hSignal(2:end,i),[timeLims yLims]);
        end
        for i=1:numFrequencyRanges
            rescalePlotHandleX(hStats(:,i),timeLims);
            rescalePlotHandleX(hSignal(1,i),timeLims);
            rescalePlotHandleX(hBursts(1,i),timeLims);


        end
    end
    function rescaleZGivenPlotHandle(plotHandles,cLims)
        [numRows,numCols] = size(plotHandles);
        for i=1:numRows
            for j=1:numCols
                clim(plotHandles(i,j),cLims);
            end
        end
    end
    function cla_Callback2(~,~)
        claGivenPlotHandle(hGridPlots2([1,2,4,5]));
    end
end

% Limits & rescaling
function rescaleGivenPlotHandle(plotHandles,axisLims)
[numRows,numCols] = size(plotHandles);
for i=1:numRows
    for j=1:numCols
        axis(plotHandles(i,j),axisLims);
    end
end
end

function rescalePlotHandleX(plotHandles,axisLims)
[numRows,numCols] = size(plotHandles);
for i=1:numRows
    for j=1:numCols
        xlim(plotHandles(i,j),axisLims);
    end
end
end


function claGivenPlotHandle(plotHandles)
[numRows,numCols] = size(plotHandles);
for i=1:numRows
    for j=1:numCols
        cla(plotHandles(i,j));
    end
end
end


