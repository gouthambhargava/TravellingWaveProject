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
burstTS = cell(numFrequencyRanges,numGoodElectrodes);
filteredSignal = cell(numFrequencyRanges,numGoodElectrodes);
if strcmp(analysisMethod,'hilbert')
    for iFreq=1:numFrequencyRanges
        for iElec=1:numGoodElectrodes
            [~,~,~,burstTS{iFreq,iElec},filteredSignal{iFreq,iElec}] = getHilbertBurst(squeeze(allData(iElec,:,:)),timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,freqRangeList{iFreq},filterOrder,1);
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
hTFAllTrials = getPlotHandles(numSelectedElectrodes,1,[0.275 0.05 0.1 0.8]);
hTFSingleTrial = getPlotHandles(numSelectedElectrodes,1,[0.4 0.05 0.15 0.8]);
hSignalSingleTrial = getPlotHandles(numSelectedElectrodes,numFrequencyRanges,[0.575 0.05 0.25 0.8]);
hBurstsAllElectrodes = getPlotHandles(numFrequencyRanges,1,[0.85 0.05 0.125 0.8],0.05,0.05);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plot_Callback(~,~)

        % Data choices
        selectedElectrodes0 = str2num(get(hElectrodes,'String')); %#ok<ST2NM>
        trialNum = trialNumList((get(hTrial,'val')));

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
        for i=1:numSelectedElectrodes
            ePos = find(selectedElectrodes(i)==goodElectrodes);
            signalAllTrials = squeeze(allData(ePos,:,:));

            if strcmp(analysisMethod,'hilbert') % Use Multi-taper method
                plotTFMT(hTFAllTrials(i),signalAllTrials,timeVals,axisRange1List,freqRangeList,colorNamesFreqRanges,'raw'); % Time-frequency power spectrum for all trials
                plotTFMT(hTFSingleTrial(i),squeeze(signalAllTrials(trialNum,:)),timeVals,axisRange1List,freqRangeList,colorNamesFreqRanges,'raw');
                hold(hTFSingleTrial(i),'on');
            else
                plotTFMP(hTFAllTrials(i),subjectName,expDate,protocolName,dataPath,gridType,sPos,oriPos,selectedElectrodes(i),timeVals,axisRange1List,freqRangeList,colorNamesFreqRanges); % Write this code - this data should be saved somewhere
            end

            % Plot filtered signal and show bursts
            for j=1:numFrequencyRanges
                plot(hSignalSingleTrial(i,j),timeVals,filteredSignal{j,ePos}(trialNum,:),'color',colorNamesFreqRanges(j,:));
                hold(hSignalSingleTrial(i,j),'on');
                plot(hSignalSingleTrial(i,j),timeVals,burstTS{j,ePos}(trialNum,:)-1,'color','r','linewidth',2);
                axis(hSignalSingleTrial(i,j),[axisRange1List{2} axisRange2List{j}]);

                plot(hTFSingleTrial(i),timeVals, mean(freqRangeList{j})+ burstTS{j,ePos}(trialNum,:)-1,'color',colorNamesFreqRanges(j,:),'linewidth',2);
            end
        end

        % plot bursts for all electrodes
        for i=1:numFrequencyRanges
            for j=1:numGoodElectrodes
                plot(hBurstsAllElectrodes(i),timeVals,j+burstTS{i,j}(trialNum,:)-1,'color','r','linewidth',1);
                hold(hBurstsAllElectrodes(i),'on');
            end
            xlim(hBurstsAllElectrodes(i),axisRange1List{2});
        end
        % % Plot pgd values with significance
        % plot(hStatsGrid(1),timeVals,outputs.pgd,'LineWidth',2)
        % xlim(hStatsGrid(1),timeRange);
    end

    function plot_Callback2(~,~)
        % Data choices
        % Filtering options
        freqRangeHz = [str2double(get(hFreqRangeMin,'String')) str2double(get(hFreqRangeMax,'String'))];
        trialNum = trialNumList((get(hTrial,'val')));
        timeRangeProp = [str2double(get(hTimeRangePropMin,'String')) str2double(get(hTimeRangePropMax,'String'))];
        timeRangeProp = dsearchn(timeVals',timeRangeProp(1)):dsearchn(timeVals',timeRangeProp(2));

        %plot options
        addPlots = get(hPlotSelect,'val');

        % Axis Ranges
        axisRange1List = cell(1,numAxisRanges1);
        for i=1:numAxisRanges1
            axisRange1List{i} = [str2double(get(hAxisRange1Min{i},'String')) str2double(get(hAxisRange1Max{i},'String'))];
        end
        axisRange2List = cell(1,numFrequencyRanges);
        for i=1:numFrequencyRanges
            axisRange2List{i} = [str2double(get(hAxisRange2Min{i},'String')) str2double(get(hAxisRange2Max{i},'String'))];
        end

        direction = outputs.direction;

        %%%%%%%%%%%%%%%%%%%%%%%%%%% Plot data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %get phi and amp grid
        [ampGrid,phiGrid,mag] = getHilbertTransData(allData(:,trialNum,:),goodElectrodes);
        [X,Y] = meshgrid(1:1:9);

        %%%%%%%%%%%% Plot single trial and indicate bursts %%%%%%%%%%%%
        xlineP1 = [];xlineP2 = [];xlineP3 = [];xlineP4 = [];xlineP5 = [];xlineP6 = [];xlineMag = [];xlinepgd = [];

        %plot phase coherence
        if addPlots==1
            plot(hStatsGrid(2),timeVals,mag,'LineWidth',1)
        else
            for gaborInd = 1:numel(goodElectrodes)
                elec = goodElectrodes(gaborInd);
                gaborPath = fullfile(dataPath,subjectName,gridType,expDate,protocolName,'mpAnalysis',[elec,num2str(selectedElectrodes(i))]);
                signalAllTrials = squeeze(allData(gaborInd,:,:));
                if ~exist(gaborPath,'file')
                    disp('Execute RunMP to run matching pursuit')
                else
                    %                       gammaAtom = zeros(1,length(goodElectrodes));
                    burstTS = nan(length(goodElectrodes),length(timeVals));
                    sigBursts = zeros(length(goodElectrodes),3);
                    load(gaborPath)
                    gaborInfo = gaborInfo(trialNum,:,:);
                    header = header(trialNum,:);

                    %initialize params for MP
                    thresholdFraction=0.5;
                    diffPower = getChangeInPower(signalAllTrials,timeVals,stimulusPeriodS,baselinePeriodS,freqRangeHz);
                    thresholdFactor = diffPower*thresholdFraction;
                    [lengthList,freqList,timeList,~,~,~] = getBurstLengthMP(signalAllTrials(trialNum,:),timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,freqRangeHz,[],[],[],gaborInfo,header);
                    if ~isempty(lengthList)
                        [~,sigBursts(gaborInd,:),burstTS(gaborInd,:), ~] = gammaBurstMPParams(lengthList,timeList,freqList,squeeze(gaborInfo),timeVals);
                    end
                end
            end
            burstTS = burstTS+(1:numel(goodElectrodes))';
            plot(hStatsGrid(2),timeVals,burstTS,'.','LineWidth',1.2,'Color','black')
            hold(hStatsGrid(2),'on')
            for ploti = 1:numel(goodElectrodes)
                plot(hStatsGrid(2),sigBursts(2),ploti,'o')
                hold(hStatsGrid(2),'on')
            end
        end

        if get(plotToggle, 'Value') == 1
            for ind = 1:numel(timeRangeProp)
                delete(xlineP1);delete(xlineP2);delete(xlineP3);delete(xlineP4);delete(xlineP5);delete(xlineP6);delete(xlineMag);delete(xlinepgd);
                xlineP1 = xline(hSignalSingleTrial(1),timeVals(timeRangeProp(ind)),'LineWidth',2);
                xlineP2 = xline(hSignalSingleTrial(2),timeVals(timeRangeProp(ind)),'LineWidth',2);
                xlineP3 = xline(hSignalSingleTrial(3),timeVals(timeRangeProp(ind)),'LineWidth',2);
                xlineP4 = xline(hSignalSingleTrial(4),timeVals(timeRangeProp(ind)),'LineWidth',2);
                xlineP5 = xline(hSignalSingleTrial(5),timeVals(timeRangeProp(ind)),'LineWidth',2);
                xlineP6 = xline(hSignalSingleTrial(6),timeVals(timeRangeProp(ind)),'LineWidth',2);
                xlineMag = xline(hStatsGrid(2), timeVals(timeRangeProp(ind)),'LineWidth',2);
                xlinepgd = xline(hStatsGrid(1), timeVals(timeRangeProp(ind)),'LineWidth',2);
                xlim(hStatsGrid(2),axisRange1List{2});

                %plot amp propagation plots
                imagesc(ampGrid(:,:,timeRangeProp(ind)),'parent',hGridPlots(3))
                caxis(hGridPlots(3),[min(ampGrid(:,:,timeRangeProp(1):timeRangeProp(end)),[],'all'),max(ampGrid(:,:,timeRangeProp(1):timeRangeProp(end)),[],'all')])
                hold(hGridPlots(3),'on')
                directions = repmat(direction(timeRangeProp(ind)),size(phiGrid,1),size(phiGrid,2));
                U = real(exp(1i*directions));
                V = imag(exp(1i*directions));
                quiver(X,Y,U,V,'Color','white','LineWidth',1,'AutoScaleFactor',0.5,'parent',hGridPlots(3))
                hold(hGridPlots(3),'off')

                %plot phase propagation plots
                imagesc(cos(phiGrid(:,:,timeRangeProp(ind))),'parent',hGridPlots(4))
                caxis(hGridPlots(4),[min(cos(phiGrid(:,:,timeRangeProp(1):timeRangeProp(end))),[],'all'),max(cos(phiGrid(:,:,timeRangeProp(1):timeRangeProp(end))),[],'all')])
                hold(hGridPlots(4),'on')
                directions = repmat(direction(timeRangeProp(ind)),size(phiGrid,1),size(phiGrid,2));
                U = real(exp(1i*directions));
                V = imag(exp(1i*directions));
                quiver(X,Y,U,V,'Color','white','LineWidth',1,'AutoScaleFactor',0.5,'parent',hGridPlots(4))
                hold(hGridPlots(4),'off')
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

        % Rescale TF plots
        rescaleGivenPlotHandle(hTFAllTrials,[timeLims freqLims]);
        rescaleGivenPlotHandle(hTFSingleTrial,[timeLims freqLims]);

        rescaleZGivenPlotHandle(hTFAllTrials,cLims);
        rescaleZGivenPlotHandle(hTFSingleTrial,cLims);

        % Rescale stats plots
        xlim(hStatsGrid(1),timeLims);
        xlim(hStatsGrid(2),timeLims);
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
        cla(hStatsGrid(2))
        cla(hGridPlots(3))
        cla(hGridPlots(4))
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

function [ampGrid,phiGrid,mag] = getHilbertTransData(signal,goodElectrodes)
%signal is a m by n matrix where m is the time points n is the number of
%channels
%check dimensions of signal and reshape if required
signal = squeeze(signal);
if size(signal,1)<size(signal,2)
    signal = signal';
end

% get phi and amp grid
gridLayout = rot90(reshape(1:81,[9,9]),2); %set grid layout alpaH
ampGrid = nan(size(gridLayout,1),size(gridLayout,2),size(signal,1));
phiGrid = nan(size(gridLayout,1),size(gridLayout,2),size(signal,1));

for i = 1:numel(goodElectrodes)
    [x,y] = find(gridLayout==goodElectrodes(i));
    ampGrid(x,y,:) = abs(hilbert(signal(:,i))).^2;
    phiGrid(x,y,:) = angle(hilbert(signal(:,i)));
end

%calculate phase coherence
%add the phases together and calculate the magnitude of this for each
%time point.
mag = zeros(1,length(phiGrid));
for ind = 1:size(phiGrid,3)
    phiValues = reshape(phiGrid(:,:,ind),[1,numel(phiGrid(:,:,ind))]);
    phiValues(isnan(phiValues)) = [];
    mag(ind) = sqrt(sum(phiValues).^2);
end
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
function plotTFMT(hTF,data,timeVals,axisRanges,freqRangeHz,colorNames,type)

Fs = round(1/(timeVals(2)-timeVals(1)));

% TF options
movingwin = [0.25 0.025];
params.tapers = [1 1];
params.pad = -1;
params.Fs = Fs;
params.fpass = [0 200];
params.trialave = 1; % averaging across trials
blRange = [-0.5 0];

%time frequency analysis with multitapers
[S,timeTF,freqVals] = mtspecgramc(data',movingwin,params);
xValToPlot = timeTF+timeVals(1)-1/Fs;
logPower = log10(S);

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
clim(hTF,axisRanges{3});

% Indicate frequency ranges of interest
numRanges = length(freqRangeHz);
for i=1:numRanges
    line([axisRanges{2}],[freqRangeHz{i}(1) freqRangeHz{i}(1)],'color',colorNames(i,:),'parent',hTF);
    line([axisRanges{2}],[freqRangeHz{i}(2) freqRangeHz{i}(2)],'color',colorNames(i,:),'parent',hTF);
end
end

function getBurstInformation
%%%%%%%%%% Parameters for burst calculation, hardcoded for now %%%%%%%%%%%%
filterOrder = 4;
thresholdFactor = 1;
baselinePeriodS = [-0.5 0];
stimulusPeriodS = [0.25 0.75];
fRange = [0 100]; cLims = [-25 25]; timeRange = [-0.5 1];
freqRangeHz = [str2double(get(hFreqRangeMin,'String')) str2double(get(hFreqRangeMax,'String'))];
trialNum = trialNumList((get(hTrial,'val')));
disp('Calculating TW parameters')
[outputs] = getTWCircParams(squeeze(allData(:,trialNum,:)),timeVals,goodElectrodes,freqRangeHz,0,[]);

end