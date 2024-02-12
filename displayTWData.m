function displayTWData3(subjectName,expDate,protocolName,dataPath,sizePos,selectedElectrodes)

fontSizeSmall = 10; fontSizeMedium = 12; fontSizeLarge = 16;
backgroundColor = 'w'; panelHeight = 0.1;
gridType = 'Microelectrode';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Getting all data...');
[allData,goodElectrodes,timeVals,rfData,parameters] = loadLFPData(subjectName,expDate,protocolName,dataPath,gridType,sizePos);

%%%%%%%%%%%%%%%%%%%%%%%% Plot RF information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hGridPlots = getPlotHandles(4,1,[0.025 0.05 0.1 0.8],0.025,0.05,0);
numSelectedElectrodes = length(selectedElectrodes);
colorNames = jet(numSelectedElectrodes);
showRFPositionsSelectedElectrodes(hGridPlots,goodElectrodes,selectedElectrodes,rfData,parameters,colorNames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Filtering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPanel1 = uipanel('Title','Filtering','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.025 1-panelHeight 0.175 panelHeight]);

% Frequency Range
freqRange0 = [30 60]; 
uicontrol('Parent',hPanel1,'Unit','Normalized','Position',[0 0.5 0.5 0.5],'Style','text','String','Freq Range (Hz)','FontSize',fontSizeMedium);
hFreqRangeMin = uicontrol('Parent',hPanel1,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.5 0.5 0.25 0.5], ...
    'Style','edit','String',num2str(freqRange0(1)),'FontSize',fontSizeSmall);
hFreqRangeMax = uicontrol('Parent',hPanel1,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.75 0.5 0.25 0.5], ...
    'Style','edit','String',num2str(freqRange0(2)),'FontSize',fontSizeSmall);

% Method
uicontrol('Parent',hPanel1,'Unit','Normalized','Position',[0 0 0.5 0.5],'Style','text','String','Filter method','FontSize',fontSizeMedium);
methodList = [{'bandpass'} {'MP'}];
hMethod = uicontrol('Parent',hPanel1,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.5 0 0.5 0.5],'Style','popup','String',methodList,'FontSize',fontSizeMedium);

%%%%%%%%%%%%%%%%%%%%%%%%%%%  Data selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPanel2 = uipanel('Title','Data','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.2 1-panelHeight 0.2 panelHeight]);

% Selected electrodes
uicontrol('Parent',hPanel2,'Unit','Normalized','Position',[0 0.5 0.25 0.5],'Style','text','String','Elecs','FontSize',fontSizeMedium);
hElectrodes = uicontrol('Parent',hPanel2,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.25 0.5 0.75 0.5],'Style','edit','String',num2str(selectedElectrodes),'FontSize',fontSizeMedium);

% Trial Number
numTrials = size(allData,2);
uicontrol('Parent',hPanel2,'Unit','Normalized','Position',[0 0 0.25 0.5],'Style','text','String','TrialNum','FontSize',fontSizeMedium);
trialNumList = 1:numTrials;
hTrial = uicontrol('Parent',hPanel2,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.25 0 0.15 0.5],'Style','popup','String',trialNumList,'FontSize',fontSizeMedium);

% Plot data
uicontrol('Parent',hPanel2,'Unit','Normalized','Position',[0.4 0 0.2 0.5],'Style','pushbutton','String','plot','FontSize',fontSizeMedium,'Callback',{@plot_Callback});
uicontrol('Parent',hPanel2,'Unit','Normalized','Position',[0.6 0 0.2 0.5],'Style','pushbutton','String','Rescale','FontSize',fontSizeMedium,'Callback',{@rescale_Callback});
uicontrol('Parent',hPanel2,'Unit','Normalized','Position',[0.8 0 0.2 0.5],'Style','pushbutton','String','Clear','FontSize',fontSizeMedium,'Callback',{@cla_Callback});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Axis Ranges %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPanel3 = uipanel('Title','AxisRanges1','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.4 1-panelHeight 0.15 panelHeight]);
axisRange1List0{1} = [0 100]; axisRange1Name{1} = 'FreqLims (Hz)';
axisRange1List0{2} = [-0.5 1]; axisRange1Name{2} = 'TimeLims (S)';

numAxisRanges1 = length(axisRange1List0);
hAxisRange1Min = cell(1,numAxisRanges1);
hAxisRange1Max = cell(1,numAxisRanges1);

for ii=1:numAxisRanges1
    uicontrol('Parent',hPanel3,'Unit','Normalized','Position',[0 1-ii/numAxisRanges1 0.5 1/numAxisRanges1],'Style','text','String',axisRange1Name{ii},'FontSize',fontSizeSmall);
    hAxisRange1Min{ii} = uicontrol('Parent',hPanel3,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.5 1-ii/numAxisRanges1 0.25 1/numAxisRanges1], ...
    'Style','edit','String',num2str(axisRange1List0{ii}(1)),'FontSize',fontSizeSmall);
    hAxisRange1Max{ii} = uicontrol('Parent',hPanel3,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.75 1-ii/numAxisRanges1 0.25 1/numAxisRanges1], ...
    'Style','edit','String',num2str(axisRange1List0{ii}(2)),'FontSize',fontSizeSmall);
end

hPanel4 = uipanel('Title','AxisRanges2','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.55 1-panelHeight 0.15 panelHeight]);
axisRange2List0{1} = [-100 100]; axisRange2Name{1} = 'YRange (S)';
axisRange2List0{2} = [-1.5 1.5]; axisRange2Name{2} = 'cLims (topo)';

numAxisRanges2 = length(axisRange2List0);
hAxisRange2Min = cell(1,numAxisRanges2);
hAxisRange2Max = cell(1,numAxisRanges2);

for ii=1:numAxisRanges2
    uicontrol('Parent',hPanel4,'Unit','Normalized','Position',[0 1-ii/numAxisRanges2 0.5 1/numAxisRanges2],'Style','text','String',axisRange2Name{ii},'FontSize',fontSizeSmall);
    hAxisRange2Min{ii} = uicontrol('Parent',hPanel4,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.5 1-ii/numAxisRanges2 0.25 1/numAxisRanges2], ...
    'Style','edit','String',num2str(axisRange2List0{ii}(1)),'FontSize',fontSizeSmall);
    hAxisRange2Max{ii} = uicontrol('Parent',hPanel4,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.75 1-ii/numAxisRanges2 0.25 1/numAxisRanges2], ...
    'Style','edit','String',num2str(axisRange2List0{ii}(2)),'FontSize',fontSizeSmall);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%  Phase propagation plot panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPanel5 = uipanel('Title','Plot Phase Propagation','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.70 1-panelHeight 0.2 panelHeight]);

timeRangeprop0 = [-0.5 1]; 
uicontrol('Parent',hPanel5,'Unit','Normalized','Position',[0 0.5 0.5 0.5],'Style','text','String','Time Range (s)','FontSize',fontSizeMedium);
hTimeRangePropMin = uicontrol('Parent',hPanel5,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.5 0.5 0.25 0.5], ...
    'Style','edit','String',num2str(timeRangeprop0(1)),'FontSize',fontSizeSmall);
hTimeRangePropMax = uicontrol('Parent',hPanel5,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.75 0.5 0.25 0.5], ...
    'Style','edit','String',num2str(timeRangeprop0(2)),'FontSize',fontSizeSmall);

plotToggle = uicontrol('Parent',hPanel5,'Unit','Normalized','Position',[0.4 0 0.4 0.5],'Style','togglebutton','String','Plot/Pause','FontSize',fontSizeMedium,'Callback',{@plot_Callback2},'Value',0);
uicontrol('Parent',hPanel5,'Unit','Normalized','Position',[0.8 0 0.2 0.5],'Style','pushbutton','String','Clear','FontSize',fontSizeMedium,'Callback',{@cla_Callback2});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot handles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hTFAllTrials = getPlotHandles(numSelectedElectrodes,1,[0.15 0.05 0.1 0.8]);
hTFSingleTrial = getPlotHandles(numSelectedElectrodes,1,[0.275 0.05 0.25 0.8]);
hSignalSingleTrial = getPlotHandles(numSelectedElectrodes,1,[0.55 0.05 0.25 0.8]);
hStatsGrid = getPlotHandles(2,1,[0.825 0.05 0.15 0.8],0.05,0.05);

%%%%%%%%%% Parameters for burst calculation, hardcoded for now %%%%%%%%%%%%
filterOrder = 4; 
thresholdFactor = 1;
baselinePeriodS = [-0.5 0]; 
stimulusPeriodS = [0.25 0.75];
fRes=2; sd=0.025; fRange = [0 100]; cLims = [-2.5 2.5]; timeRange = [-0.5 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plot_Callback(~,~)

        % Filtering options
        freqRangeHz = [str2double(get(hFreqRangeMin,'String')) str2double(get(hFreqRangeMax,'String'))];
        filteringMethod = get(hMethod,'val');

        % Data choices
        selectedElectrodes0 = str2num(get(hElectrodes,'String')); %#ok<ST2NM>
        trialNum = trialNumList((get(hTrial,'val')));

        if ~isequal(selectedElectrodes0,selectedElectrodes)
            showRFPositionsSelectedElectrodes(hGridPlots,goodElectrodes,selectedElectrodes0,rfData,parameters,colorNames);
            selectedElectrodes = selectedElectrodes0;
        end

        % Axis Ranges
        axisRange1List = cell(1,numAxisRanges1);
        for i=1:numAxisRanges1
            axisRange1List{i} = [str2double(get(hAxisRange1Min{i},'String')) str2double(get(hAxisRange1Max{i},'String'))];
        end
        axisRange2List = cell(1,numAxisRanges2);
        for i=1:numAxisRanges2
            axisRange2List{i} = [str2double(get(hAxisRange2Min{i},'String')) str2double(get(hAxisRange2Max{i},'String'))];
        end

        % calculate TW params for plotting
        [outputs] = getTWCircParams(squeeze(allData(:,trialNum,:)),timeVals,goodElectrodes,[30 60],0);
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%% Plot data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i=1:numSelectedElectrodes
            signalAllTrials = squeeze(allData(selectedElectrodes(i),:,:));
            cgt = zeros(numel(fRange(1):fRes:fRange(2)),size(signalAllTrials,2),size(signalAllTrials,1));
            
              
            if filteringMethod==1 % Bandpass filter
                [~,burstStartS,burstEndS,burstTS,bpfSignalAllTrials,hilbertPowerAllTrials] = getBurstLengthHilbert(signalAllTrials(trialNum,:),timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,freqRangeHz,filterOrder,1);
                for tfi = 1:size(signalAllTrials,1)
                    [cgt(:,:,tfi),freqVals]=getCGT(signalAllTrials(tfi,:),timeVals,fRange,fRes,sd);
                end
                TFPow = log10(abs(cgt).^2);   
                basePowS = repmat(mean(TFPow(:,(dsearchn(timeVals',baselinePeriodS(1)):dsearchn(timeVals',baselinePeriodS(2))),trialNum),2),1,size(TFPow,2));
                basePowA = repmat(mean(mean(TFPow(:,(dsearchn(timeVals',baselinePeriodS(1)):dsearchn(timeVals',baselinePeriodS(2))),:),2),3),1,size(TFPow,2));
                corrPowS = TFPow(:,:,trialNum)- basePowS;
                corrPowA = mean(TFPow,3) - basePowA;

             elseif filteringMethod==2 %matching pursuit
                 [~,~,~,gaborInfo,header,~] = getBurstLengthMP(signalAllTrials,timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,freqRangeHz);   
                 
                 %reconstruct the signal and energy from the gamma atoms
                 scale=1;freq=2;pos=3;amp=5;phi=6;
                 for triali = 1:size(gaborInfo,1)   
                    [energyT(:,:,triali),freqVals]=getEnergyMP3p1(gaborInfo(triali,:,:),header(triali,:),timeVals);
                 end
                 
                 TFRange = find(freqVals<100);
                 freqVals = freqVals(TFRange);
                 corrPowS = energyT(TFRange,:,trialNum);
                 corrPowA = mean(energyT(TFRange,:,:),3);
                 gaborInfo = squeeze(gaborInfo(trialNum,:,:));
                 reconMP = [];
                 header = header(trialNum,:);
                 for gAtoms = 1:size(gaborInfo,1) 
                     reconMP(gAtoms,:) = gabor(header(2)/header(1),header(1),gaborInfo(gAtoms,scale),gaborInfo(gAtoms,freq),gaborInfo(gAtoms,pos),gaborInfo(gAtoms,amp),gaborInfo(gAtoms,phi));
                 end
                 reqFreqs = find(gaborInfo(:,2)>30 & gaborInfo(:,2)<60);   
                 reconMP = sum(reconMP(reqFreqs,:));
                 
                 % Get burst info via hilbert 
                 [~,burstStartS,burstEndS,burstTS,bpfSignalAllTrials,hilbertPowerAllTrials] = getBurstLengthHilbert(reconMP,timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,freqRangeHz,filterOrder,0);        
            end

            %%%%%%%%%%%% Plot single trial and indicate bursts %%%%%%%%%%%%  
            %plot the burst
            %Uncomment to plot the power time series 
%             plot(hSignalSingleTrial(i),timeVals,zscore(hilbertPowerAllTrials),'color',colorNames(i,:),'linewidth',0.5);
%             hold(hSignalSingleTrial(i),'on'); 
            plot(hSignalSingleTrial(i),timeVals,bpfSignalAllTrials,'color',colorNames(i,:),'linewidth',2);
            hold(hSignalSingleTrial(i),'on');
            plot(hSignalSingleTrial(i),timeVals,burstTS,'color','k','linewidth',2);
            hold(hSignalSingleTrial(i),'off');
            
            % plot TF single trial
            pcolor(hTFSingleTrial(i),timeVals,freqVals,corrPowS);
            shading(hTFSingleTrial(i),'interp');
            colormap(hTFSingleTrial(i), 'jet')
            axis(hTFSingleTrial(i),[timeRange fRange]);
            caxis(hTFSingleTrial(i),cLims);
            
            % plot TF all trial
            pcolor(hTFAllTrials(i),timeVals,freqVals,corrPowA);
            shading(hTFAllTrials(i),'interp');
            colormap(hTFAllTrials(i), 'jet')        
            axis(hTFAllTrials(i),[timeRange fRange]);
            caxis(hTFAllTrials(i),cLims);
        end
         % Plot pgd values with significance
         plot(hStatsGrid(1),timeVals,outputs.pgd,'LineWidth',2)
         hold(hStatsGrid(1),'on')
         yline(hStatsGrid(1),0.5,'LineWidth',0.75,'Color','red')
         hold(hStatsGrid(1),'off')
         xlim(hStatsGrid(1),timeRange);
        
    end
    function cla_Callback(~,~)
        claGivenPlotHandle(hTFAllTrials);
        claGivenPlotHandle(hTFSingleTrial);
        claGivenPlotHandle(hSignalSingleTrial);
        claGivenPlotHandle(hStatsGrid);
        cla(hStatsGrid(1))
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

    function plot_Callback2(~,~)
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % Filtering options
         freqRangeHz = [str2double(get(hFreqRangeMin,'String')) str2double(get(hFreqRangeMax,'String'))];
         % Data choices
         selectedElectrodes0 = str2num(get(hElectrodes,'String')); %#ok<ST2NM>
         trialNum = trialNumList((get(hTrial,'val')));
         timeRangeProp = [str2double(get(hTimeRangePropMin,'String')) str2double(get(hTimeRangePropMax,'String'))];
         timeRangeProp = dsearchn(timeVals',timeRangeProp(1)):dsearchn(timeVals',timeRangeProp(2));
         filteringMethod = get(hMethod,'val');
         
         % Axis Ranges
         axisRange1List = cell(1,numAxisRanges1);
         for i=1:numAxisRanges1
            axisRange1List{i} = [str2double(get(hAxisRange1Min{i},'String')) str2double(get(hAxisRange1Max{i},'String'))];
         end
         axisRange2List = cell(1,numAxisRanges2);
         for i=1:numAxisRanges2
            axisRange2List{i} = [str2double(get(hAxisRange2Min{i},'String')) str2double(get(hAxisRange2Max{i},'String'))];
         end
         
         %Calculate direction params
         [outputs] = getTWCircParams(squeeze(allData(:,trialNum,:)),timeVals,goodElectrodes,[30 60],0);
         direction = outputs.direction;
         %%%%%%%%%%%%%%%%%%%%%%%%%%% Plot data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if filteringMethod==1
             req = 1;
         else 
             req = 2;
         end
         for i=1:numSelectedElectrodes
            signalAllTrials = squeeze(allData(selectedElectrodes(i),:,:));
            [~,burstStartS{i},burstEndS{i},burstTS{i},bpfSignalAllTrials{i},~] = getBurstLengthHilbert(signalAllTrials,timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,freqRangeHz,filterOrder,req);          
         end 
         
                  
         %get phi and amp grid
         [ampGrid,phiGrid,mag] = getHilbertTransData(allData(:,trialNum,:),timeVals,freqRangeHz,4,goodElectrodes);
         [X,Y] = meshgrid(1:1:9);
         
         %%%%%%%%%%%% Plot single trial and indicate bursts %%%%%%%%%%%%
         xlineP1 = [];xlineP2 = [];xlineP3 = [];xlineP4 = [];xlineP5 = [];xlineP6 = [];xlineMag = [];xlinepgd = [];
             
         %plot phase coherence
         
         plot(hStatsGrid(2), timeVals,mag,'LineWidth',1)
         
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
                    U = cos(phiGrid(:,:,timeRangeProp(ind))); 
                    V = sin(phiGrid(:,:,timeRangeProp(ind))); 
                    caxis(hGridPlots(3),[min(ampGrid(:,:,timeRangeProp(1):timeRangeProp(end)),[],'all'),max(ampGrid(:,:,timeRangeProp(1):timeRangeProp(end)),[],'all')])
                    hold(hGridPlots(3),'on')
                    quiver(X,Y,U,V,'Color','white','LineWidth',1,'AutoScaleFactor',0.5,'parent',hGridPlots(3))
                
                    %plot phase propagation plots
                    imagesc(cos(phiGrid(:,:,timeRangeProp(ind))),'parent',hGridPlots(4))
                    caxis(hGridPlots(4),[min(cos(phiGrid(:,:,timeRangeProp(1):timeRangeProp(end))),[],'all'),max(cos(phiGrid(:,:,timeRangeProp(1):timeRangeProp(end))),[],'all')])
                    hold(hGridPlots(4),'on')
                    U1 = direction(timeRangeProp(ind))*X;
                    V1 = direction(timeRangeProp(ind))*Y;
                    quiver(X,Y,U1,V1,'parent',hGridPlots(4))
                    hold(hGridPlots(4),'off')
                    drawnow

                end 
                else
                uiwait
            end
    end

    function cla_Callback2(~,~)
        
        cla(hStatsGrid(2))
        cla(hGridPlots(3))
        cla(hGridPlots(4))
    end
end


function [allData,goodElectrodes,timeVals,rfData,parameters] = loadLFPData(subjectName,expDate,protocolName,dataPath,gridType,sizePos)

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
goodPos = setdiff(parameters.parameterCombinations{1,1,sizePos,1,size(parameters.parameterCombinations,5)},badTrials);

allData = zeros(numGoodElectrodes,length(goodPos),length(timeVals));

for i=1:numGoodElectrodes
    lfpData = load(fullfile(folderName,'segmentedData','LFP',['elec' num2str(goodElectrodes(i)) '.mat']));
    allData(i,:,:) = lfpData.analogData(goodPos,:);
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

