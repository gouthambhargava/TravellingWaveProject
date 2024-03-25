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

function displayTWData(subjectName,expDate,protocolName,dataPath,sPos,oriPos,selectedElectrodes,analysisMethod,freqRangeList)

fontSizeSmall = 10; fontSizeMedium = 12; fontSizeLarge = 16;
backgroundColor = 'w'; panelHeight = 0.125;
gridType = 'Microelectrode';

numFrequencyRanges = length(freqRangeList);
colorNamesFreqRanges = gray(numFrequencyRanges+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Getting data...');
[allData,goodElectrodes,timeVals,rfData,parameters] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sPos,oriPos);
numGoodElectrodes = length(goodElectrodes);


%%%%%%%%%%%%%%%%%%%%%%% Do burst estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thresholdFactor = 1;
baselinePeriodS = [-0.5 0];
stimulusPeriodS = [0.25 0.75];
filterOrder = 4;

disp('Performing burst analysis...');
burstTS = zeros(numGoodElectrodes,size(allData,2),size(allData,3),numFrequencyRanges);
filteredSignal = zeros(numGoodElectrodes,size(allData,2),size(allData,3),numFrequencyRanges);
if strcmp(analysisMethod,'hilbert')
    for iFreq=1:numFrequencyRanges
        for iElec=1:numGoodElectrodes
            [~,~,~,burstTS(iElec,:,:,iFreq),filteredSignal(iElec,:,:,iFreq)] = getHilbertBurst(squeeze(allData(iElec,:,:)),timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,freqRangeList{iFreq},filterOrder,1);
        end
    end
else
end

%%%%%%%%%%%%%%%%%%%%%%%% Plot RF information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hGridPlots = getPlotHandles(3,2,[0.025 0.05 0.225 0.8],0.025,0.05,0);
numSelectedElectrodes = length(selectedElectrodes);
colorNamesElectrodes = jet(numSelectedElectrodes);
showRFPositionsSelectedElectrodes(hGridPlots(1,:),goodElectrodes,selectedElectrodes,rfData,parameters,colorNamesElectrodes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%  Data selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPanel1 = uipanel('Title','Data','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0 1-panelHeight 0.2 panelHeight]);

% Selected electrodes
uicontrol('Parent',hPanel1,'Unit','Normalized','Position',[0 0.5 0.25 0.5],'Style','text','String','Elecs','FontSize',fontSizeMedium);
hElectrodes = uicontrol('Parent',hPanel1,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.25 0.5 0.75 0.5],'Style','edit','String',num2str(selectedElectrodes),'FontSize',fontSizeMedium);

% Trial Number
numTrials = size(allData,2);
uicontrol('Parent',hPanel1,'Unit','Normalized','Position',[0 0 0.25 0.5],'Style','text','String','TrialNum','FontSize',fontSizeMedium);
trialNumList = 1:numTrials;
hTrial = uicontrol('Parent',hPanel1,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.25 0 0.15 0.5],'Style','popup','String',trialNumList,'FontSize',fontSizeMedium);

% Plot data
uicontrol('Parent',hPanel1,'Unit','Normalized','Position',[0.4 0 0.2 0.5],'Style','pushbutton','String','plot','FontSize',fontSizeMedium,'Callback',{@plot_Callback});
uicontrol('Parent',hPanel1,'Unit','Normalized','Position',[0.6 0 0.2 0.5],'Style','pushbutton','String','Rescale','FontSize',fontSizeMedium,'Callback',{@rescale_Callback});
uicontrol('Parent',hPanel1,'Unit','Normalized','Position',[0.8 0 0.2 0.5],'Style','pushbutton','String','Clear','FontSize',fontSizeMedium,'Callback',{@cla_Callback});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Axis Ranges %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPanel2 = uipanel('Title','AxisRanges1','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.2 1-panelHeight 0.2 panelHeight]);
axisRange1List0{1} = [0 100]; axisRange1Name{1} = 'FreqLims (Hz)';
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

hPanel3 = uipanel('Title','AxisRanges2','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.4 1-panelHeight 0.2 panelHeight]);

axisRange2List0 = cell(1,numFrequencyRanges);
axisRange2Name = cell(1,numFrequencyRanges);
for iFreq=1:numFrequencyRanges
    axisRange2List0{iFreq} = [-100 100]; 
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

%%%%%%%%%%%%%%%%%%%%  Phase propagation plot panel %%%%%%%%%%%%%%%%%%%%%%%%
hPanel4 = uipanel('Title','Plot Phase Propagation','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.6 1-panelHeight 0.2 panelHeight]);

timeRangeprop0 = [-0.5 1];
uicontrol('Parent',hPanel4,'Unit','Normalized','Position',[0 0.5 0.5 0.5],'Style','text','String','Time Range (s)','FontSize',fontSizeMedium);
hTimeRangePropMin = uicontrol('Parent',hPanel4,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.5 0.5 0.25 0.5], ...
    'Style','edit','String',num2str(timeRangeprop0(1)),'FontSize',fontSizeSmall);
hTimeRangePropMax = uicontrol('Parent',hPanel4,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.75 0.5 0.25 0.5], ...
    'Style','edit','String',num2str(timeRangeprop0(2)),'FontSize',fontSizeSmall);

plotToggle = uicontrol('Parent',hPanel4,'Unit','Normalized','Position',[0 0 0.5 0.5],'Style','togglebutton','String','Plot/Pause','FontSize',fontSizeMedium,'Callback',{@plot_Callback2},'Value',0);
uicontrol('Parent',hPanel4,'Unit','Normalized','Position',[0.5 0 0.5 0.5],'Style','pushbutton','String','Clear','FontSize',fontSizeMedium,'Callback',{@cla_Callback2});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot handles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hSignalSingleTrial = getPlotHandles(numSelectedElectrodes+1,numFrequencyRanges,[0.575 0.05 0.25 0.8]); %use the top most panels to plot both phsae coh and pgd
hBurstsAllElectrodes = getPlotHandles(numFrequencyRanges,1,[0.85 0.05 0.125 0.65],0.05,0.05);
hPGDPanel = getPlotHandles(2,1,[0.85 0.75 0.125 0.25],0.01,0.01);
hTFAllTrials = getPlotHandles(numSelectedElectrodes+1,1,[0.275 0.05 0.1 0.8]);
hTFSingleTrial = getPlotHandles(numSelectedElectrodes+1,1,[0.4 0.05 0.15 0.8]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get TW params %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Getting TW parameters')
trialNum = trialNumList((get(hTrial,'val')));
[outputsSG] = getTWCircParams(squeeze(allData(:,trialNum,:)),timeVals,goodElectrodes,freqRangeList{1},1,[]);
[outputsFG] = getTWCircParams(squeeze(allData(:,trialNum,:)),timeVals,goodElectrodes,freqRangeList{2},1,[]);
magSG = outputsSG.mag;
magFG = outputsFG.mag;
durIndicesSG = nan(size(magSG));
durIndicesSG(islocalmin(magSG,'MinSeparation',50)) = magSG(islocalmin(magSG,'MinSeparation',50));

durIndicesFG = nan(size(magFG));
durIndicesFG(islocalmin(magFG,'MinSeparation',25)) = magFG(islocalmin(magFG,'MinSeparation',25));
        
pgdSG = outputsSG.pgd;
pgdFG = outputsFG.pgd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function plot_Callback(~,~)

        % Data choices
        selectedElectrodes0 = str2num(get(hElectrodes,'String')); %#ok<ST2NM>
        
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
        plotTFMT(hTFAllTrials(1),allData,timeVals,axisRange1List,freqRangeList,colorNamesFreqRanges,'raw',1,[],[]); % Time-frequency power spectrum for all trials
        hold(hTFAllTrials(1),'on')
        xline(hTFAllTrials(1),0.25);
        hold(hTFAllTrials(1),'on');
        xline(hTFAllTrials(1),0.75);
        plotTFMT(hTFSingleTrial(1),allData,timeVals,axisRange1List,freqRangeList,colorNamesFreqRanges,'raw',0,trialNum,[]);
        hold(hTFSingleTrial(1),'on')
        xline(hTFSingleTrial(1),0.25);
        hold(hTFSingleTrial(1),'on');
        xline(hTFSingleTrial(1),0.75);
        
        for i=1:numSelectedElectrodes
            ePos = find(selectedElectrodes(i)==goodElectrodes);
            signalAllTrials = allData(ePos,:,:);

            if strcmp(analysisMethod,'hilbert') % Use Multi-taper method
                plotTFMT(hTFAllTrials(i+1),signalAllTrials,timeVals,axisRange1List,freqRangeList,colorNamesFreqRanges,'raw',1,[],[]); % Time-frequency power spectrum for all trials
                plotTFMT(hTFSingleTrial(i+1),signalAllTrials,timeVals,axisRange1List,freqRangeList,colorNamesFreqRanges,'raw',0,trialNum,[]);
                hold(hTFSingleTrial(i+1),'on');
%                 hold(hTFAllTrials(i+1),'on');
%                 xline(hTFAllTrials(i+1),0.25);
%                 hold(hTFAllTrials(i+1),'on');
%                 xline(hTFAllTrials(i+1),0.75);
                
            else
                plotTFMP(hTFAllTrials(i),subjectName,expDate,protocolName,dataPath,gridType,sPos,oriPos,selectedElectrodes(i),timeVals,axisRange1List,freqRangeList,colorNamesFreqRanges); % Write this code - this data should be saved somewhere
            end

            % Plot filtered signal and show bursts
            for j=1:numFrequencyRanges
                plot(hSignalSingleTrial(i+1,j),timeVals,squeeze(filteredSignal(ePos,trialNum,:,j)),'color',colorNamesFreqRanges(j,:));
                hold(hSignalSingleTrial(i+1,j),'on');
                plot(hSignalSingleTrial(i+1,j),timeVals,squeeze(burstTS(ePos,trialNum,:,j))-1,'color','r','linewidth',2);
                axis(hSignalSingleTrial(i+1,j),[axisRange1List{2} axisRange2List{j}]);
                plot(hTFSingleTrial(i+1),timeVals, mean(freqRangeList{j})+ squeeze(burstTS(ePos,trialNum,:,j))-1,'color',colorNamesFreqRanges(j,:),'linewidth',2);
%                 hold(hSignalSingleTrial(i+1,j),'on');
%                 xline(hTFSingleTrial(i+1),0.25);
%                 hold(hTFSingleTrial(i+1,j),'on');
%                 xline(hTFSingleTrial(i+1),0.75);
            end
        end
       % plot bursts 
       chanVals = repmat((1:numel(goodElectrodes))',1,size(burstTS,3))';  
       burstDurationSG = timeVals(durIndicesSG>1);
       burstDurationFG = timeVals(durIndicesFG>1);
       plot(timeVals,squeeze(burstTS(:,trialNum,:,1))'+chanVals-1,'parent',hBurstsAllElectrodes(1),'color','r','linewidth',1);
       hold(hBurstsAllElectrodes(1),'on')
       
       for xlineInd = 1:length(burstDurationSG)   
           xline(burstDurationSG(xlineInd),'parent',hBurstsAllElectrodes(1),'color','k');
           hold(hBurstsAllElectrodes(1),'on')
       end
       
       plot(timeVals,squeeze(burstTS(:,trialNum,:,2))'+chanVals-1,'parent',hBurstsAllElectrodes(2),'color','r','linewidth',1);
       hold(hBurstsAllElectrodes(2),'on')   
       for xlineInd = 1:length(burstDurationFG)   
           xline(burstDurationFG(xlineInd),'parent',hBurstsAllElectrodes(2),'color','k');
           hold(hBurstsAllElectrodes(2),'on')
       end
       
        %plot phase coherence
         plot(hSignalSingleTrial(1,1),timeVals,magSG,'color',colorNamesFreqRanges(1,:))
         axis(hSignalSingleTrial(1,1),[axisRange1List{2} axisRange2List{j}]);
%        xlim(hSignalSingleTrial(1,1),[timeVals(1) timeVals(end)]);
         hold(hSignalSingleTrial(1,1),'on')
         plot(hSignalSingleTrial(1,1),timeVals,durIndicesSG,'*','Color','red')
         hold(hSignalSingleTrial(1,1),'off')
            
         plot(hSignalSingleTrial(1,2),timeVals,magFG,'color',colorNamesFreqRanges(2,:))
         hold(hSignalSingleTrial(1,2),'on')
         plot(hSignalSingleTrial(1,2),timeVals,durIndicesFG,'*','Color','red')
         axis(hSignalSingleTrial(1,2),[axisRange1List{2} axisRange2List{j}]);
%        xlim(hSignalSingleTrial(1,2),[timeVals(1) timeVals(end)]);
         hold(hSignalSingleTrial(1,2),'off')
         
        %plot pgd
        plot(hPGDPanel(1),timeVals,pgdSG,'color',colorNamesFreqRanges(1,:))
        hold(hPGDPanel(1),'on')
        plot(hPGDPanel(2),timeVals,pgdFG,'color',colorNamesFreqRanges(2,:))
        hold(hPGDPanel(2),'on')
        xlim(hPGDPanel(1),axisRange1List{2});
        ylim(hPGDPanel(1),[-1 1])
        xlim(hPGDPanel(2),axisRange1List{2});
        ylim(hPGDPanel(2),[-1 1])
%         xlim(hPGDPanel(1),[timeVals(1) timeVals(end)]);
%         xlim(hPGDPanel(2),[timeVals(1) timeVals(end)]);
    end

    function plot_Callback2(~,~)
        % Data choices
        % Filtering options
        trialNum = trialNumList((get(hTrial,'val')));
        timeRangeProp = [str2double(get(hTimeRangePropMin,'String')) str2double(get(hTimeRangePropMax,'String'))];
        timeRangeProp = dsearchn(timeVals',timeRangeProp(1)):dsearchn(timeVals',timeRangeProp(2));

        % Axis Ranges
        axisRange1List = cell(1,numAxisRanges1);
        for i=1:numAxisRanges1
            axisRange1List{i} = [str2double(get(hAxisRange1Min{i},'String')) str2double(get(hAxisRange1Max{i},'String'))];
        end
        axisRange2List = cell(1,numFrequencyRanges);
        for i=1:numFrequencyRanges
            axisRange2List{i} = [str2double(get(hAxisRange2Min{i},'String')) str2double(get(hAxisRange2Max{i},'String'))];
        end

        directionSG = outputsSG.direction;
        directionFG = outputsFG.direction;

        %%%%%%%%%%%%%%%%%%%%%%%%%%% Plot data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %get phi and amp grid
        
        [ampGridSG,phiGridSG] = getHilbertTransData(filteredSignal(:,trialNum,:,1),goodElectrodes);
        [ampGridFG,phiGridFG] = getHilbertTransData(filteredSignal(:,trialNum,:,2),goodElectrodes);
        [X,Y] = meshgrid(1:1:9);

        %%%%%%%%%%%% Plot single trial and indicate bursts %%%%%%%%%%%%
        xlineP1 = [];xlineP2 = [];xlineP3 = [];xlineP4 = [];xlineP5 = [];xlineP6 = [];xlineMag1 = [];xlineMag2 = [];xlinePGD1 = [];xlinePGD2 = [];xlineBurst1 = [];xlineBurst2 = [];
    
        if get(plotToggle, 'Value') == 1
            for ind = 1:numel(timeRangeProp)
                delete(xlineP1);delete(xlineP2);delete(xlineP3);delete(xlineP4);delete(xlineP5);delete(xlineP6);delete(xlineMag1);delete(xlineMag2);delete(xlinePGD1);delete(xlinePGD2);delete(xlineBurst1);delete(xlineBurst2);
                xlineP1 = xline(hSignalSingleTrial(1),timeVals(timeRangeProp(ind)),'LineWidth',2);
                xlineP2 = xline(hSignalSingleTrial(2),timeVals(timeRangeProp(ind)),'LineWidth',2);
                xlineP3 = xline(hSignalSingleTrial(3),timeVals(timeRangeProp(ind)),'LineWidth',2);
                xlineP4 = xline(hSignalSingleTrial(4),timeVals(timeRangeProp(ind)),'LineWidth',2);
                xlineP5 = xline(hSignalSingleTrial(5),timeVals(timeRangeProp(ind)),'LineWidth',2);
                xlineP6 = xline(hSignalSingleTrial(6),timeVals(timeRangeProp(ind)),'LineWidth',2);
                xlineMag1 = xline(hSignalSingleTrial(1,1), timeVals(timeRangeProp(ind)),'LineWidth',2);
                xlineMag2 = xline(hSignalSingleTrial(1,2), timeVals(timeRangeProp(ind)),'LineWidth',2);
                xlinePGD1 = xline(hPGDPanel(1), timeVals(timeRangeProp(ind)),'LineWidth',2);
                xlinePGD2 = xline(hPGDPanel(2), timeVals(timeRangeProp(ind)),'LineWidth',2);
                xlineBurst1 = xline(hBurstsAllElectrodes(1), timeVals(timeRangeProp(ind)),'LineWidth',2);
                xlineBurst2 = xline(hBurstsAllElectrodes(2), timeVals(timeRangeProp(ind)),'LineWidth',2);
                
                
                %plot amp propagation plots
                % For slow gamma
                imagesc(ampGridSG(:,:,timeRangeProp(ind)),'parent',hGridPlots(2))
                caxis(hGridPlots(2),[min(ampGridSG(:,:,timeRangeProp(1):timeRangeProp(end)),[],'all'),max(ampGridSG(:,:,timeRangeProp(1):timeRangeProp(end)),[],'all')])
                hold(hGridPlots(2),'on')
                directions = repmat(directionSG(timeRangeProp(ind)),size(phiGridSG,1),size(phiGridSG,2));
                U = real(exp(1i*directions));
                V = imag(exp(1i*directions));
                quiver(X,Y,U,V,'Color','white','LineWidth',1,'AutoScaleFactor',0.5,'parent',hGridPlots(2))
                hold(hGridPlots(2),'off')
                 
                % For fast gamma
                imagesc(ampGridFG(:,:,timeRangeProp(ind)),'parent',hGridPlots(5))
                caxis(hGridPlots(5),[min(ampGridFG(:,:,timeRangeProp(1):timeRangeProp(end)),[],'all'),max(ampGridFG(:,:,timeRangeProp(1):timeRangeProp(end)),[],'all')])
                hold(hGridPlots(5),'on')
                directions = repmat(directionFG(timeRangeProp(ind)),size(phiGridFG,1),size(phiGridFG,2));
                U = real(exp(1i*directions));
                V = imag(exp(1i*directions));
                 quiver(X,Y,U,V,'Color','white','LineWidth',1,'AutoScaleFactor',0.5,'parent',hGridPlots(5))
                 hold(hGridPlots(5),'off')
                
                %plot phase propagation plots
                %for slow gamma
                imagesc(cos(phiGridSG(:,:,timeRangeProp(ind))),'parent',hGridPlots(3,1))
                caxis(hGridPlots(3,1),[min(cos(phiGridSG(:,:,timeRangeProp(1):timeRangeProp(end))),[],'all'),max(cos(phiGridSG(:,:,timeRangeProp(1):timeRangeProp(end))),[],'all')])
                hold(hGridPlots(3,1),'on')
                directions = repmat(directionSG(timeRangeProp(ind)),size(phiGridSG,1),size(phiGridSG,2));
                U = real(exp(1i*directions));
                V = imag(exp(1i*directions));
                quiver(X,Y,U,V,'Color','white','LineWidth',1,'AutoScaleFactor',0.5,'parent',hGridPlots(3,1))
                hold(hGridPlots(3,1),'off')
                
                %for fast gamma
                imagesc(cos(phiGridFG(:,:,timeRangeProp(ind))),'parent',hGridPlots(3,2))
                caxis(hGridPlots(3,2),[min(cos(phiGridFG(:,:,timeRangeProp(1):timeRangeProp(end))),[],'all'),max(cos(phiGridFG(:,:,timeRangeProp(1):timeRangeProp(end))),[],'all')])
                hold(hGridPlots(3,2),'on')
                directions = repmat(directionFG(timeRangeProp(ind)),size(phiGridFG,1),size(phiGridFG,2));
                U = real(exp(1i*directions));
                V = imag(exp(1i*directions));
                quiver(X,Y,U,V,'Color','white','LineWidth',1,'AutoScaleFactor',0.5,'parent',hGridPlots(3,2))
                hold(hGridPlots(3,2),'off')
                
                drawnow
            end
        else
            uiwait
        end
    end


    function cla_Callback(~,~)
        claGivenPlotHandle(hTFAllTrials);
        claGivenPlotHandle(hTFSingleTrial);
        claGivenPlotHandle(hSignalSingleTrial);
        claGivenPlotHandle(hBurstsAllElectrodes);
        claGivenPlotHandle(hPGDPanel)
        

        function claGivenPlotHandle(plotHandles)
            [numRows,numCols] = size(plotHandles);
            for i=1:numRows
                for j=1:numCols
                    cla(plotHandles(i,j));
                end
            end
        end
    end

    function rescale_Callback(~,~)
        freqLims = [str2double(get(hAxisRange1Min{1},'String')) str2double(get(hAxisRange1Max{1},'String'))];
        timeLims = [str2double(get(hAxisRange1Min{2},'String')) str2double(get(hAxisRange1Max{2},'String'))];
        yLims = [str2double(get(hAxisRange2Min{1},'String')) str2double(get(hAxisRange2Max{1},'String'))];
        cLims = [str2double(get(hAxisRange2Min{2},'String')) str2double(get(hAxisRange2Max{2},'String'))];

        rescaleGivenPlotHandle(hSignalSingleTrial,[timeLims yLims]);
        
        %rescale pgd plots 
        xlim(hPGDPanel(1),timeLims);
        xlim(hPGDPanel(2),timeLims);
        
        %rescale burst plots
        xlim(hBurstsAllElectrodes(1),timeLims);
        xlim(hBurstsAllElectrodes(2),timeLims);
        
        % Rescale TF plots
        rescaleGivenPlotHandle(hTFAllTrials,[timeLims freqLims]);
        rescaleGivenPlotHandle(hTFSingleTrial,[timeLims freqLims]);

%         rescaleZGivenPlotHandle(hTFAllTrials,cLims);
%         rescaleZGivenPlotHandle(hTFSingleTrial,cLims);

%         % Rescale stats plots
%         xlim(hStatsGrid(1),timeLims);
%         xlim(hStatsGrid(2),timeLims);
    end

    function rescaleZGivenPlotHandle(plotHandles,cLims)
        [numRows,numCols] = size(plotHandles);
        for i=1:numRows
            for j=1:numCols
                caxis(plotHandles(i,j),cLims);
            end
        end
    end



    function cla_Callback2(~,~)
        cla(hGridPlots(2));
        cla(hGridPlots(3));
        cla(hGridPlots(5));
        cla(hGridPlots(6));
    end
end

function stimElectrode = showElectrodeRFs(hPlot,goodElectrodes,colorName,rfData,parameters)
aValsUnique = parameters.aValsUnique;
eValsUnique = parameters.eValsUnique;

N = length(goodElectrodes);
rfAziList = zeros(1,N);
rfEleList = zeros(1,N);

axes(hPlot);
hold(hPlot,'on');
for i=1:N
    rfAziList(i) = rfData.rfStats(goodElectrodes(i)).meanAzi;
    rfEleList(i) = rfData.rfStats(goodElectrodes(i)).meanEle;
    plot(hPlot,rfAziList(i),rfEleList(i),'marker','+','color',colorName);
    text(rfAziList(i),rfEleList(i),num2str(goodElectrodes(i)),'color',colorName,'fontsize',8);
end
plot(aValsUnique,eValsUnique,'rs','markerfacecolor','r');

distanceList = sqrt((aValsUnique-rfAziList).^2+(eValsUnique-rfEleList).^2);
stimElectrode = goodElectrodes(distanceList==min(distanceList));
end
function electrodeArray = showElectrodePositions(hPlot,highlightElectrodes,colorNames,hideElectrodeNums)

if ~exist('highlightElectrodes','var'); highlightElectrodes=[];         end
if ~exist('colorNames','var');          colorNames=[];                  end
if ~exist('hideElectrodeNums','var');    hideElectrodeNums=0;           end

electrodeArray = ...
    [81 72 63 54 45 36 27 18 09;
    80 71 62 53 44 35 26 17 08;
    79 70 61 52 43 34 25 16 07;
    78 69 60 51 42 33 24 15 06;
    77 68 59 50 41 32 23 14 05;
    76 67 58 49 40 31 22 13 04;
    75 66 57 48 39 30 21 12 03;
    74 65 56 47 38 29 20 11 02;
    73 64 55 46 37 28 19 10 01];

[numRows,numCols] = size(electrodeArray);

axes(hPlot);
dX = 1/numCols;
dY = 1/numRows;

lineXRow = zeros(2,numRows);lineYRow = zeros(2,numRows);
for i=1:numRows
    lineXRow(:,i) = [0 1]; lineYRow(:,i) = [i*dY i*dY];
end
lineXCol = zeros(2,numCols);lineYCol = zeros(2,numCols);
for i=1:numCols
    lineXCol(:,i) = [i*dX i*dX]; lineYCol(:,i) = [0 1];
end
line(lineXRow,lineYRow,'color','k'); hold on;
line(lineXCol,lineYCol,'color','k');
hold off;

if ~isempty(highlightElectrodes)
    for i=1:length(highlightElectrodes)
        highlightElectrode=highlightElectrodes(i);

        [highlightRow,highlightCol] = find(highlightElectrode==electrodeArray);

        % Create patch
        patchX = (highlightCol-1)*dX;
        patchY = (numRows-highlightRow)*dY;
        patchLocX = [patchX patchX patchX+dX patchX+dX];
        patchLocY = [patchY patchY+dY patchY+dY patchY];

        if iscell(colorNames)
            patch('XData',patchLocX,'YData',patchLocY,'FaceColor',colorNames{i});
        else
            patch('XData',patchLocX,'YData',patchLocY,'FaceColor',colorNames);
        end
    end
end

if ~hideElectrodeNums
    % Write electrode numbers
    for i=1:numRows
        textY = (numRows-i)*dY + dY/2;
        for j=1:numCols
            textX = (j-1)*dX + dX/2;
            if electrodeArray(i,j)>0
                text(textX,textY,num2str(electrodeArray(i,j)),'HorizontalAlignment','center');
            end
        end
    end
end

set(hPlot,'XTickLabel',[],'YTickLabel',[]);
end

function showRFPositionsSelectedElectrodes(hRFPlots,goodElectrodes,selectedElectrodes,rfData,parameters,colorNames)
cla(hRFPlots(1)); cla(hRFPlots(2));
showElectrodeRFs(hRFPlots(1),goodElectrodes,'k',rfData,parameters);

electrodeArray = showElectrodePositions(hRFPlots(2),[],[]);
showElectrodePositions(hRFPlots(2),setdiff(electrodeArray,goodElectrodes),'k');

for i=1:length(selectedElectrodes)
    showElectrodeRFs(hRFPlots(1),selectedElectrodes(i),colorNames(i,:),rfData,parameters);
    showElectrodePositions(hRFPlots(2),selectedElectrodes(i),colorNames(i,:));
end
xlabel(hRFPlots(1),'Azimuth (deg)');
ylabel(hRFPlots(1),'Elevation (deg)');
end

function [ampGrid,phiGrid] = getHilbertTransData(signal,goodElectrodes)
%signal is a m by n matrix where m is the time points n is the number of
%channels
%check dimensions of signal and reshape if required
signal = squeeze(signal);

% get phi and amp grid
gridLayout = rot90(reshape(1:81,[9,9]),2); %set grid layout alpaH
ampGrid = nan(size(gridLayout,1),size(gridLayout,2),size(signal,2));
phiGrid = nan(size(gridLayout,1),size(gridLayout,2),size(signal,2));

for i = 1:numel(goodElectrodes)
    [x,y] = find(gridLayout==goodElectrodes(i));
    ampGrid(x,y,:) = abs(hilbert(signal(i,:))).^2;
    phiGrid(x,y,:) = angle(hilbert(signal(i,:)));
end

%calculate phase coherence
%add the phases together and calculate the magnitude of this for each
%time point.
% mag = zeros(1,length(phiGrid));
% for ind = 1:size(phiGrid,3)
%     phiValues = reshape(phiGrid(:,:,ind),[1,numel(phiGrid(:,:,ind))]);
%     phiValues(isnan(phiValues)) = [];
%     mag(ind) = sqrt(sum(phiValues.^2));
% end
clear x y
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
function yLims = getYLims(plotHandles)

[numRows,numCols] = size(plotHandles);
% Initialize
yMin = inf;
yMax = -inf;

for row=1:numRows
    for column=1:numCols
        % get positions
        axis(plotHandles(row,column),'tight');
        tmpAxisVals = axis(plotHandles(row,column));
        if tmpAxisVals(3) < yMin
            yMin = tmpAxisVals(3);
        end
        if tmpAxisVals(4) > yMax
            yMax = tmpAxisVals(4);
        end
    end
end

yLims=[yMin yMax];
end

% Time-frequency analysis
function plotTFMT(hTF,data,timeVals,axisRanges,freqRangeHz,colorNames,type,trialAvg,trialNo,chan)
% supply data in the form channels x trails x times
% specify if trial avegaring is required, 0 not required, 1 required
if nargin<10
    chan=[];
end
Fs = round(1/(timeVals(2)-timeVals(1)));

if trialAvg==0
    data = data(:,trialNo,:);
end

% TF options
movingwin = [0.25 0.025];
params.tapers = [1 1];
params.pad = -1;
params.Fs = Fs;
params.fpass = [0 200];
params.trialave = trialAvg; % averaging across trials
blRange = [-0.5 0];


for i = 1:size(data,1)
%time frequency analysis with multitapers
[S,timeTF,freqVals] = mtspecgramc(squeeze(data(i,:,:))',movingwin,params);
xValToPlot = timeTF+timeVals(1)-1/Fs;
logPower(:,:,i) = log10(S);
end

if ~isempty(chan)
    logPower = logPower(:,:,chan);
else 
    logPower = mean(logPower,3);
end

if strcmp(type,'delta')
    blPos = intersect(find(xValToPlot>=blRange(1)),find(xValToPlot<blRange(2)));
    logBLPower = repmat(mean(logPower(blPos,:,:),1),length(xValToPlot),1);
    deltaPower = 10*(logPower - logBLPower);
    pcolor(hTF,xValToPlot,freqVals,deltaPower');
else
    pcolor(hTF,xValToPlot,freqVals,logPower');
end

shading(hTF,'interp');
colormap(hTF,'jet');
axis(hTF,[axisRanges{2} axisRanges{1}]);
caxis(hTF,axisRanges{3});
stimValues = [dsearchn(xValToPlot',0.25),dsearchn(xValToPlot',0.75)];

% Indicate frequency ranges of interest
numRanges = length(freqRangeHz);
for i=1:numRanges
    line([axisRanges{2}],[freqRangeHz{i}(1) freqRangeHz{i}(1)],'color',colorNames(i,:),'parent',hTF);
    line([axisRanges{2}],[freqRangeHz{i}(2) freqRangeHz{i}(2)],'color',colorNames(i,:),'parent',hTF);
end

end