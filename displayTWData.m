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
thresholdFactor = 4;
stimulusDurationS = [0 0.8]; % Stimulus duration to be highlighted
baselinePeriodS = [-0.5 0];
stimulusPeriodS = [0.25 0.75];
analysisPeriodS = [-0.5 1];
filterOrder = 4;

disp('Performing burst analysis...');
burstTS = zeros(numGoodElectrodes,size(allData,2),size(allData,3),numFrequencyRanges);
filteredSignal = zeros(numGoodElectrodes,size(allData,2),size(allData,3),numFrequencyRanges);
if strcmp(analysisMethod,'hilbert')
    for iFreq=1:numFrequencyRanges
        for iElec=1:numGoodElectrodes
            [~,~,~,burstTS(iElec,:,:,iFreq),filteredSignal(iElec,:,:,iFreq)] = getHilbertBurst(squeeze(allData(iElec,:,:)),timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,freqRangeList{iFreq},filterOrder,1,analysisPeriodS);
        end
    end
else
    %add matching pursuit 
end

%%%%%%%%%%%%%%%%%%%%%%%% Plot RF information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hGridPlots = getPlotHandles(1,2,[0.025 0.65 0.3 0.2],0.01,0.02,0);
numSelectedElectrodes = length(selectedElectrodes);
colorNamesElectrodes = jet(numSelectedElectrodes);
electrodeArray = showRFPositionsSelectedElectrodes(hGridPlots,goodElectrodes,selectedElectrodes,rfData,parameters,colorNamesElectrodes);

% get location for each electrode for segmentation into clusters
locList = nan(numGoodElectrodes,2);
for iElec = 1:numGoodElectrodes
    [locList(iElec,1),locList(iElec,2)] = find(electrodeArray==goodElectrodes(iElec));
end

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

hPanel3 = uipanel('Title','AxisRanges2','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.4 1-panelHeight 0.2 panelHeight]);

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
%this needs to be still be implemented
hPanel4 = uipanel('Title','Traveling Wave','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.6 1-panelHeight 0.2 panelHeight]);
uicontrol('Parent',hPanel4,'Unit','Normalized','Position',[0 0.5 0.25 0.5],'Style','text','String','Electrode Fraction','FontSize',fontSizeSmall);
hElecFrac = uicontrol('Parent',hPanel4,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.25 0.7 0.2 0.3], ...
    'Style','edit','String','0.6','FontSize',fontSizeSmall);

uicontrol('Parent',hPanel4,'Unit','Normalized','Position',[0 0 0.25 0.5],'Style','text','String','Select Electrodes','FontSize',fontSizeSmall);
hSelectElec = uicontrol('Parent',hPanel4,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.25 0.2 0.2 0.3], ...
    'Style','popup','String',{'All','Selected'},'FontSize',fontSizeSmall);

uicontrol('Parent',hPanel4,'Unit','Normalized','Position',[0.5 0.5 0.25 0.5],'Style','text','String','Ref Choice','FontSize',fontSizeSmall);
hSelectRef = uicontrol('Parent',hPanel4,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.75 0.72 0.2 0.3], ...
    'Style','popup','String',{'Avg',num2str(goodElectrodes')},'FontSize',fontSizeSmall);

% electrodeFraction = 0.6;
% electrodeChoice = 'all'; % 'all' or 'selected'
% refPhaseChoice = 'avg'; % choose electrode number

%%%%%%%%%%%%%%%%%%%%  Phase propagation plot panel %%%%%%%%%%%%%%%%%%%%%%%%
hPanel5 = uipanel('Title','Plot Phase Propagation','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.8 1-panelHeight 0.2 panelHeight]);

timeRangeprop0 = [-0.5 1];
uicontrol('Parent',hPanel5,'Unit','Normalized','Position',[0 0.5 0.5 0.5],'Style','text','String','Time Range (s)','FontSize',fontSizeMedium);
hTimeRangePropMin = uicontrol('Parent',hPanel5,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.5 0.5 0.25 0.5], ...
    'Style','edit','String',num2str(timeRangeprop0(1)),'FontSize',fontSizeSmall);
hTimeRangePropMax = uicontrol('Parent',hPanel5,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.75 0.5 0.25 0.5], ...
    'Style','edit','String',num2str(timeRangeprop0(2)),'FontSize',fontSizeSmall);

plotToggle = uicontrol('Parent',hPanel5,'Unit','Normalized','Position',[0 0 0.5 0.5],'Style','togglebutton','String','Plot/Pause','FontSize',fontSizeMedium,'Callback',{@plot_Callback2},'Value',0);
uicontrol('Parent',hPanel5,'Unit','Normalized','Position',[0.5 0 0.5 0.5],'Style','pushbutton','String','Clear','FontSize',fontSizeMedium,'Callback',{@cla_Callback2});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot handles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hTF = getPlotHandles(numSelectedElectrodes+1,2,[0.025 0.05 0.3 0.55],0.025,0.025);
hSignal = getPlotHandles(numSelectedElectrodes+1,numFrequencyRanges,[0.35 0.05 0.425 0.55],0.025,0.025);
hStats = getPlotHandles(3,numFrequencyRanges,[0.35 0.625 0.425 0.225],0.025,0.01);
hGridPlots2 = getPlotHandles(3,numFrequencyRanges,[0.8 0.05 0.19 0.8]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Needs to be plotted only once
plotTFMT(hTF(1,1),allData,timeVals,axisRange1List0,freqRangeList,colorNamesFreqRanges,'raw',[],stimulusDurationS); % Time-frequency power spectrum for all trials

phaseMatrix = cell(1,numFrequencyRanges);
outputsTW = cell(1,numFrequencyRanges);
state = 0; % state variable to display the phases
directions = cell(1,numFrequencyRanges);

    function plot_Callback(~,~)

        % Data choices
        selectedElectrodes0 = str2num(get(hElectrodes,'String')); %#ok<ST2NM>
        trialNum = trialNumList((get(hTrial,'val')));
        
        %get burst fraction statistics
        elecFrac = str2double(get(hElecFrac,'string'));
        if get(hSelectElec,'val') == 2
            elecFrac = 0;
        end
        
        for i = 1:numFrequencyRanges
            burstFrac{i} = nansum(squeeze(burstTS(:,trialNum,:,i)));
            sigValues = burstFrac{i}>numGoodElectrodes*elecFrac;
            burstFrac{i}(sigValues) = 1;
            burstFrac{i}(setdiff(1:length(burstFrac{i}),find(sigValues))) = nan;
        end

        % Calculate TW parameters
        disp('Getting TW parameters');

        phaseMatrix = cell(1,numFrequencyRanges);
        burstMatrix = cell(1,numFrequencyRanges);
        outputsTW = cell(1,numFrequencyRanges);

        for i = 1:numFrequencyRanges
            phaseMatrix{i} = angle(hilbert(squeeze(filteredSignal(:,trialNum,:,i))'))';
            tmp = squeeze(burstTS(:,trialNum,:,i));
            tmp(isnan(tmp)) = 0;
            burstMatrix{i} = tmp;
            outputsTW{i} = getTWCircParams(phaseMatrix{i},burstMatrix{i},timeVals,goodElectrodes,locList,elecFrac,'all');
            directions{i} = outputsTW{i}.direction;
        end

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
                plot(hSignal(i+1,j),timeVals,squeeze(filteredSignal(ePos,trialNum,:,j)),'color',colorNamesFreqRanges(j,:));
                hold(hSignal(i+1,j),'on');
                plot(hSignal(i+1,j),timeVals,squeeze(burstTS(ePos,trialNum,:,j))-1,'color','r','linewidth',2);
                axis(hSignal(i+1,j),[axisRange1List{2} axisRange2List{j}]);
                plot(hTF(i+1,2),timeVals, mean(freqRangeList{j})+ squeeze(burstTS(ePos,trialNum,:,j))-1,'color',colorNamesFreqRanges(j,:),'linewidth',2);
            end
        end

            

        % Plot TW stats
        for i = 1:numFrequencyRanges
            
            %get burst durations
            burstSeg = burstFrac{i};
            burstSeg(isnan(burstSeg)) = 0;
            indices = find(burstSeg==0);
            dirLimits1 = find(diff(indices)>1)';
            dirLimits2 = [indices(dirLimits1)+1;indices(dirLimits1+1)-1]';
            dirLimits2(find((dirLimits2(:,2)-dirLimits2(:,1))<10),:) = [];
            colors = parula(size(dirLimits2,1));
            
            %Plot bursts
            chanVals = repmat((1:numel(goodElectrodes))',1,size(burstTS,3))';
            plot(timeVals,squeeze(burstTS(:,trialNum,:,i))'+chanVals-1,'parent',hSignal(1,i),'color','r','linewidth',1);
            hold(hSignal(1,i),'on');
            for j = 1:size(dirLimits2,1)
                burstTemp = nan(1,length(timeVals));
                burstTemp(dirLimits2(j,1):dirLimits2(j,2)) = 1;
                plot(timeVals,burstTemp*79,'|','parent',hSignal(1,i),'color',colors(j,:),'linewidth',1);
                hold(hSignal(1,i),'on');
            end
            xlim(hSignal(1,i),axisRange1List{2});
            
            plot(hStats(1,i),timeVals,angle(exp(1i*outputsTW{i}.direction)),'color',colorNamesFreqRanges(i,:));
            axis(hStats(1,i),[axisRange1List{2} -pi pi]);

            plot(hStats(2,i),timeVals,outputsTW{i}.Wavelength,'color',colorNamesFreqRanges(i,:));
            xlim(hStats(2,i),axisRange1List{2});

            plot(hStats(3,i),timeVals,outputsTW{i}.coh,'color','r');
            hold(hStats(3,i),'on');
            plot(hStats(3,i),timeVals,abs(mean(exp(1i*phaseMatrix{i}))),'color','m');
            plot(hStats(3,i),timeVals,mean(burstMatrix{i}),'color','k');
            plot(hStats(3,i),timeVals,outputsTW{i}.pgd,'color','g');
            hold(hStats(3,i),'on');
%           yline(hStats(3,i),elecFrac,'color','b');
            line(axisRange1List{2},[elecFrac elecFrac],'parent',hStats(3,i),'color','b');
            axis(hStats(3,i),[axisRange1List{2} 0 1]);
            
            %Plot angle plots for calculated directions
            for duri = 1:size(dirLimits2,1)
                polarPlotTW(hGridPlots2(3,i),directions{i}(dirLimits2(duri,1):dirLimits2(duri,2)),10,colors(duri,:));
            end   
        end
    end

    function plot_Callback2(~,~)
        
           % get inputs
           refPhaseChoice = get(hSelectRef,'val'); 
           [X,Y] = meshgrid(1:1:9);

        if get(plotToggle,'Value') == 1

            state = 1; % This allows the program to start displaying the phase plots

            timeRangeProp = [str2double(get(hTimeRangePropMin,'String')) str2double(get(hTimeRangePropMax,'String'))];
            timeRangeProp = dsearchn(timeVals',timeRangeProp(1)):dsearchn(timeVals',timeRangeProp(2));

            %%%%%%%%%%%% Plot single trial and indicate bursts %%%%%%%%%%%%
            % Initialize
            for j=1:numFrequencyRanges
                for k=1:size(hSignal,1)
                    xLinesSignal(k,j) = line([timeVals(timeRangeProp(1)) timeVals(timeRangeProp(1))],get(hSignal(k,j),'YLim'),'parent',hSignal(k,j),'LineWidth',2); %#ok<*AGROW>
                end
                for k=1:size(hStats,1)
                    xLinesStats(k,j) = line([timeVals(timeRangeProp(1)) timeVals(timeRangeProp(1))],get(hStats(k,j),'YLim'),'parent',hStats(k,j),'LineWidth',2);
                end
            end

            for ind = 1:numel(timeRangeProp)
                if state==0 % if state changes to 0, delete lines and return
                    for j=1:numFrequencyRanges
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

                    % Delete old lines and create new ones
                    for k=1:size(hSignal,1)
                        delete(xLinesSignal(k,j));
                        xLinesSignal(k,j) = line([timeVals(timeRangeProp(ind)) timeVals(timeRangeProp(ind))],get(hSignal(k,j),'YLim'),'parent',hSignal(k,j),'LineWidth',2);
                    end
                    for k=1:size(hStats,1)
                        delete(xLinesStats(k,j));
                        xLinesStats(k,j) = line([timeVals(timeRangeProp(ind)) timeVals(timeRangeProp(ind))],get(hStats(k,j),'YLim'),'parent',hStats(k,j),'LineWidth',2);
                    end

                    % plot absolute phases
                    tmpPhases = phaseMatrix{j}(:,timeRangeProp(ind));
                    if refPhaseChoice==1
                        refPhase = circ_mean(tmpPhases);
                    else
                        refPhase = tmpPhases(refPhaseChoice-1,:);
                    end

                    for e=1:numGoodElectrodes
                        tmpMatrix1(locList(e,1),locList(e,2)) = tmpPhases(e);
                        tmpMatrix2(locList(e,1),locList(e,2)) = angle(exp(1i*(tmpPhases(e) - refPhase)));
                    end

                    imagesc(tmpMatrix1,'parent',hGridPlots2(1,j));
                    colormap(hGridPlots2(1,j),'parula');
                    caxis(hGridPlots2(1,j),[-pi pi]);
                    colorbar(hGridPlots2(1,j),'northoutside');

                    imagesc(tmpMatrix2,'parent',hGridPlots2(2,j));
                    colormap(hGridPlots2(2,j),'hsv');
                    caxis(hGridPlots2(2,j),[-pi pi]);
                    colorbar(hGridPlots2(2,j),'northoutside');
                    hold(hGridPlots2(2,j),'on')
                    directionGrid = repmat(directions{j}(timeRangeProp(ind)),size(tmpMatrix2,1),size(tmpMatrix2,2));
                    U = real(exp(1i*directionGrid));
                    V = imag(exp(1i*directionGrid));
                    quiver(X,Y,U,V,'Color','white','LineWidth',1,'AutoScaleFactor',0.5,'parent',hGridPlots2(2,j))
%                   hold(hGridPlots(2,1),'off')
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
        claGivenPlotHandle(hGridPlots2)
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
        claGivenPlotHandle(hGridPlots2(1:4));
    end
end

% Electrode location and display
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

function electrodeArray = showRFPositionsSelectedElectrodes(hRFPlots,goodElectrodes,selectedElectrodes,rfData,parameters,colorNames)
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

% Limits & rescaling
function rescaleGivenPlotHandle(plotHandles,axisLims)
[numRows,numCols] = size(plotHandles);
for i=1:numRows
    for j=1:numCols
        axis(plotHandles(i,j),axisLims);
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

% Time-frequency analysis
function plotTFMT(hTF,data,timeVals,axisRanges,freqRangeHz,colorNames,type,trialNo,stimulusPeriodS)
% supply data in the form channels x trials x times
% set trialNo to [] if all channels are needed

Fs = round(1/(timeVals(2)-timeVals(1)));

if ~isempty(trialNo)
    data = data(:,trialNo,:);
end

% TF options
movingwin = [0.25 0.025];
params.tapers = [1 1];
params.pad = -1;
params.Fs = Fs;
params.fpass = [0 200];
params.trialave = 1; % averaging across trials
blRange = [-0.5 0];

% Time frequency analysis with multitapers
for j = 1:size(data,1)
    [S,timeTF,freqVals] = mtspecgramc(squeeze(data(j,:,:))',movingwin,params);
    logPower(:,:,j) = log10(S);
end
xValToPlot = timeTF+timeVals(1)-1/Fs;

logPower = mean(logPower,3);

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

% Indicate frequency ranges of interest
numRanges = length(freqRangeHz);
for j=1:numRanges
    line([axisRanges{2}],[freqRangeHz{j}(1) freqRangeHz{j}(1)],'color',colorNames(j,:),'parent',hTF);
    line([axisRanges{2}],[freqRangeHz{j}(2) freqRangeHz{j}(2)],'color',colorNames(j,:),'parent',hTF);
end

if ~isempty(stimulusPeriodS)
    line([stimulusPeriodS(1) stimulusPeriodS(1)],[axisRanges{1}],'color','k','parent',hTF);
    line([stimulusPeriodS(2) stimulusPeriodS(2)],[axisRanges{1}],'color','k','parent',hTF);
end
end

function polarPlotTW(axes,phaseValues,binWidth,color)

%bin the data and get counts
nBins = 360/binWidth;
allBins = linspace(0,360,nBins);
wrappedAngles = wrapTo360(rad2deg(phaseValues));
binnedVals = discretize(wrappedAngles,allBins);
binnedVals(isnan(binnedVals)) = [];
uniqueVals = unique(binnedVals);

counts = zeros(1,numel(uniqueVals));
for i = 1:numel(uniqueVals)
    counts(uniqueVals(i)) = numel(find(binnedVals==uniqueVals(i)));
end
counts = counts/max(counts); %normalize counts 

%plot the counts
%get bin co-ordinates
binCord = exp(1i*deg2rad([allBins,allBins(2)]));

for i = 1:numel(counts)
    xcords = counts(i)*[real(binCord(i)),real(binCord(i+1))];
    ycords = counts(i)*[imag(binCord(i)),imag(binCord(i+1))];
    patch([0 xcords 0],[0 ycords 0],color,'parent',axes)
    hold(axes,'on')
end

line([-1 1],[0 0],'Color',[0.811764705882353,0.811764705882353,0.811764705882353],'parent',axes)
line([0 0],[-1 1],'Color',[0.811764705882353,0.811764705882353,0.811764705882353],'parent',axes)
circle(axes,0,0,1,[0,0,0]);
circle(axes,0,0,0.5,[0.811764705882353,0.811764705882353,0.811764705882353]);
axis(axes,[-1,1,-1,1])
text(1.01,0,['0',char(176)],'parent',axes)
text(0,1.05,['90',char(176)],'parent',axes)
text(-1.2,0,['180',char(176)],'parent',axes)
text(0,-1.05,['270',char(176)],'parent',axes)
axis(axes,'square')
axis(axes,'off')

    function circle(axes,x,y,r,colorVal)
        hold(axes,'on')
        th = 0:pi/50:2*pi;
        xunit = r * cos(th) + x;
        yunit = r * sin(th) + y;
        plot(axes,xunit, yunit,'Color',colorVal);
        hold(axes,'off')
    end
end
