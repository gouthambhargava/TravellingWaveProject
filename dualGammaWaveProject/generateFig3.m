%% calculate and plot wave directions - overlapping and nonoverlapping
% get overlapping waves
% outputsM1 is a cell of TW data for all trials for both frequency
% range (slow and fast gamma) in different rows for subject 1. Outputs for ori60deg and
% sf of 1 cpm are being considered here. 
numTrials = numel(outputsM1(1,:)); % same number of trials in both output structures

% define some parameters for wave detection
lengthLimit = 25; %ms
boundryLims = [0.25 0.75];
wobbleLim = 5; %degree
segOption = 3;
numFrequencyRanges = size(outputsM1,1);

%initialize outputs
waveVector1 = nan(numTrials,length(timeVals),numFrequencyRanges);
uniqueDirs1 = cell(numFrequencyRanges,numTrials);
waveBounds1 = cell(numFrequencyRanges,numTrials);
for i = 1:numTrials
    for j = 1:numFrequencyRanges
    [waveVector1(i,:,j),uniqueDirs1{j,i},waveBounds1{j,i}] = getWaveSegments(outputsM1{j,i},timeVals,wobbleLim,segOption,boundryLims, lengthLimit);
    end
end

% find overlapping waves 
allUniqueDirs1 = [];
dirSG1 = nan(numTrials,length(timeVals));
dirFG1 = nan(numTrials,length(timeVals));
overlap = 0.5;
for i = 1:numTrials
    [~,dirSG1(i,:),dirFG1(i,:),uniqueDirs] = getOverlappingWaves(waveVector1(i,:,1),waveBounds1{1,i},waveVector1(i,:,2),waveBounds1{2,i},overlap);
    allUniqueDirs1 = cat(2,allUniqueDirs1,uniqueDirs);
end
% get circular correlation for M1
[rho1, pval1] = circ_corrcc(allUniqueDirs1(1,:),allUniqueDirs1(2,:));

%% calculate the same parameters for M2
numTrials = numel(outputsM2(1,:)); % same number of trials in both output structures
%initialize outputs
waveVector2 = nan(numTrials,length(timeVals),numFrequencyRanges);
uniqueDirs2 = cell(numFrequencyRanges,numTrials);
waveBounds2 = cell(numFrequencyRanges,numTrials);
for i = 1:numTrials
    for j = 1:numFrequencyRanges
    [waveVector2(i,:,j),uniqueDirs2{j,i},waveBounds2{j,i}] = getWaveSegments(outputsM2{j,i},timeVals,wobbleLim,segOption,boundryLims, lengthLimit);
    end
end

% find overlapping waves 
allUniqueDirs2 = [];
dirSG2 = nan(numTrials,length(timeVals));
dirFG2 = nan(numTrials,length(timeVals));
overlap = 0.5;
for i = 1:numTrials
    [~,dirSG2(i,:),dirFG2(i,:),uniqueDirs] = getOverlappingWaves(waveVector2(i,:,1),waveBounds2{1,i},waveVector2(i,:,2),waveBounds2{2,i},overlap);
    allUniqueDirs2 = cat(2,allUniqueDirs2,uniqueDirs);
end
% get circular correlation for M1
[rho2, pval2] = circ_corrcc(allUniqueDirs2(1,:),allUniqueDirs2(2,:));

%% plot direction
binWidth = 10;
tiledlayout(2,3,'TileSpacing','Compact');
colorVals = cat(1,[52 148 186]./255,[236 112 22]./255);
h1 = nexttile;
makePolarPlot({cell2mat(uniqueDirs1(1,:)),cell2mat(uniqueDirs1(2,:))},binWidth,h1,colorVals)
title('Directions - Unsorted')

h2 = nexttile;
makePolarPlot({allUniqueDirs1(1,:),allUniqueDirs1(2,:)},binWidth,h2,colorVals)
title('Directions - Overlapping')

% get scatter plot
nexttile
x = cell2mat(dirsOL1(1,:));
y = cell2mat(dirsOL1(2,:));
scatter(x,y,'filled','LineWidth',0.01)
xlabel('Slow Gamma')
ylabel('Fast Gamma')
title('SG/FG Directions')
xticks([-3.14 0 3.14])
xticklabels({'-180','0','180'})
yticks([-3.14 0 3.14])
yticklabels({'-180','0','180'})
annotation('textbox', [0.78, 0.8, 0.1, 0.1], 'String', ['r:',num2str(round(abs(rho1),4)),', p:',num2str(round(pval1,3))])
xlim([-3.14 3.14])
ylim([-3.14 3.14])
axis square

h3 = nexttile;
makePolarPlot({cell2mat(uniqueDirs2(1,:)),cell2mat(uniqueDirs2(2,:))},binWidth,h3,colorVals)
title('Directions - Unsorted')
h4 = nexttile;
makePolarPlot({allUniqueDirs2(1,:),allUniqueDirs2(2,:)},binWidth,h4,colorVals)
title('Directions - Overlapping')

% get scatter plot
nexttile
x = cell2mat(dirsOL2(1,:));
y = cell2mat(dirsOL2(2,:));
scatter(x,y,'filled','LineWidth',0.01)
xlabel('Slow Gamma')
ylabel('Fast Gamma')
title('SG/FG Directions')
xlim([-3.14 3.14])
ylim([-3.14 3.14])
xticks([-3.14 0 3.14])
xticklabels({'-180','0','180'})
yticks([-3.14 0 3.14])
yticklabels({'-180','0','180'})
annotation('textbox', [0.71, 0.08, 0.1, 0.1], 'String', ['r:',num2str(round(rho2,3)),', p:',num2str(round(pval2,3))])
axis square
%%
annotation('textbox', [0.04, 0.565, 0.2, 0.2], 'String', 'M1','FontSize',24,'FontWeight','bold','EdgeColor','none')
annotation('textbox', [0.04, 0.11, 0.2, 0.2], 'String', 'M2','FontSize',24,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.025, 0.95, 0.02, 0.02], 'string', 'A','FontSize',24,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.3, 0.965, 0.02, 0.02], 'string', 'B','FontSize',24,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.65, 0.965, 0.02, 0.02], 'string', 'C','FontSize',24,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.025, 0.5, 0.02, 0.02], 'string', 'D','FontSize',24,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.3, 0.5, 0.02, 0.02], 'string', 'E','FontSize',24,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.65, 0.5, 0.02, 0.02], 'string', 'F','FontSize',24,'FontWeight','bold','EdgeColor','none')
