%% calculate and plot wave directions - overlapping and nonoverlapping
load('alpaH_42_0.5T_selected_met1.mat')
outputsM1 = outputs;
load('kesariH_42_0.5T_selected_met1.mat')
outputsM2 = outputs;

% define some parameters for wave detection
wobbleLim1 = 0; %degree
segOption1 = 3;
wobbleLim2 = 5; %degree
segOption2 = 2;
wobbleLim3 = 10; %degree
lengthLimit = 10; %ms
boundryLims = [0.25 0.75];
numFrequencyRanges = size(outputsM1,1);

% For M1
numTrials = numel(outputsM1(1,:)); % same number of trials in both output structures

%for 0 deg
%initialize outputs
waveVector = nan(numTrials,length(timeVals),numFrequencyRanges);
uniqueDirs = cell(numFrequencyRanges,numTrials);
waveBounds = cell(numFrequencyRanges,numTrials);

for i = 1:numTrials
    for j = 1:numFrequencyRanges
    [waveVector(i,:,j),uniqueDirs{j,i},waveBounds{j,i}] = getWaveSegments(outputsM1{j,i},timeVals,wobbleLim1,segOption1,boundryLims, lengthLimit);
    end
end

% find overlapping waves 
allUniqueDirs = [];
dirSG = nan(numTrials,length(timeVals));
dirFG = nan(numTrials,length(timeVals));
waveBoundsOv = cell(1,numTrials);
emptyCells = zeros(1,numTrials);
overlap = 0.5;
for i = 1:numTrials
    [waveBoundsOv{1,i},dirSG(i,:),dirFG(i,:),uniqueDirsTemp,emptyCells(i)] = getOverlappingWaves(waveVector(i,:,1),waveBounds{1,i},waveVector(i,:,2),waveBounds{2,i},overlap);   
    allUniqueDirs = cat(2,allUniqueDirs,uniqueDirsTemp);
end
allUniqueDirs(:,isnan(allUniqueDirs(1,:))) = [];

% get all wave angles
waveBoundsOv(emptyCells==1) = [];
ovTrials = find(emptyCells==0);

allSGWaves = [];
allFGWaves = [];
for i = 1:length(waveBoundsOv)
    waveTemp = waveBoundsOv{i};
    for j = 1:size(waveTemp{1},2)
        overlappingPts = intersect(waveTemp{1}(1,j):waveTemp{1}(2,j),waveTemp{2}(1,j):waveTemp{2}(2,j));
        sgWaves = {dirSG(ovTrials(i),overlappingPts)};
        allSGWaves = cat(1,allSGWaves,sgWaves);
        fgWaves = {dirFG(ovTrials(i),overlappingPts)};
        allFGWaves = cat(1,allFGWaves,fgWaves);
    end
end 
allOvWaves = cat(2,allSGWaves,allFGWaves);

m1_0deg = {uniqueDirs,allUniqueDirs,allOvWaves,ovTrials};    
% get circular correlation for M1
[rho1(1), pval1(1)] = circ_corrcc(allUniqueDirs(1,:),allUniqueDirs(2,:));

% for 5 deg
%initialize outputs
waveVector = nan(numTrials,length(timeVals),numFrequencyRanges);
uniqueDirs = cell(numFrequencyRanges,numTrials);
waveBounds = cell(numFrequencyRanges,numTrials);

for i = 1:numTrials
    for j = 1:numFrequencyRanges
    [waveVector(i,:,j),uniqueDirs{j,i},waveBounds{j,i}] = getWaveSegments(outputsM1{j,i},timeVals,wobbleLim2,segOption2,boundryLims, lengthLimit);
    end
end

% find overlapping waves 
allUniqueDirs = [];
dirSG = nan(numTrials,length(timeVals));
dirFG = nan(numTrials,length(timeVals));
waveBoundsOv = cell(1,numTrials);
emptyCells = zeros(1,numTrials);
overlap = 0.5;
for i = 1:numTrials
    [waveBoundsOv{1,i},dirSG(i,:),dirFG(i,:),uniqueDirsTemp,emptyCells(i)] = getOverlappingWaves(waveVector(i,:,1),waveBounds{1,i},waveVector(i,:,2),waveBounds{2,i},overlap);   
    allUniqueDirs = cat(2,allUniqueDirs,uniqueDirsTemp);
end
allUniqueDirs(:,isnan(allUniqueDirs(1,:))) = [];

% get all wave angles
waveBoundsOv(emptyCells==1) = [];
ovTrials = find(emptyCells==0);

allSGWaves = [];
allFGWaves = [];
for i = 1:length(waveBoundsOv)
    waveTemp = waveBoundsOv{i};
    for j = 1:size(waveTemp{1},2)
        overlappingPts = intersect(waveTemp{1}(1,j):waveTemp{1}(2,j),waveTemp{2}(1,j):waveTemp{2}(2,j));
        sgWaves = {dirSG(ovTrials(i),overlappingPts)};
        allSGWaves = cat(1,allSGWaves,sgWaves);
        fgWaves = {dirFG(ovTrials(i),overlappingPts)};
        allFGWaves = cat(1,allFGWaves,fgWaves);
    end
end 
allOvWaves = cat(2,allSGWaves,allFGWaves);

m1_5deg = {uniqueDirs,allUniqueDirs,allOvWaves,ovTrials};    
% get circular correlation for M1
[rho1(2), pval1(2)] = circ_corrcc(allUniqueDirs(1,:),allUniqueDirs(2,:));

% 10 degrees
%initialize outputs
waveVector = nan(numTrials,length(timeVals),numFrequencyRanges);
uniqueDirs = cell(numFrequencyRanges,numTrials);
waveBounds = cell(numFrequencyRanges,numTrials);

for i = 1:numTrials
    for j = 1:numFrequencyRanges
    [waveVector(i,:,j),uniqueDirs{j,i},waveBounds{j,i}] = getWaveSegments(outputsM1{j,i},timeVals,wobbleLim3,segOption2,boundryLims, lengthLimit);
    end
end

% find overlapping waves 
allUniqueDirs = [];
dirSG = nan(numTrials,length(timeVals));
dirFG = nan(numTrials,length(timeVals));
waveBoundsOv = cell(1,numTrials);
emptyCells = zeros(1,numTrials);
overlap = 0.5;
for i = 1:numTrials
    [waveBoundsOv{1,i},dirSG(i,:),dirFG(i,:),uniqueDirsTemp,emptyCells(i)] = getOverlappingWaves(waveVector(i,:,1),waveBounds{1,i},waveVector(i,:,2),waveBounds{2,i},overlap);   
    allUniqueDirs = cat(2,allUniqueDirs,uniqueDirsTemp);
end
allUniqueDirs(:,isnan(allUniqueDirs(1,:))) = [];

% get all wave angles
waveBoundsOv(emptyCells==1) = [];
ovTrials = find(emptyCells==0);

allSGWaves = [];
allFGWaves = [];
for i = 1:length(waveBoundsOv)
    waveTemp = waveBoundsOv{i};
    for j = 1:size(waveTemp{1},2)
        overlappingPts = intersect(waveTemp{1}(1,j):waveTemp{1}(2,j),waveTemp{2}(1,j):waveTemp{2}(2,j));
        sgWaves = {dirSG(ovTrials(i),overlappingPts)};
        allSGWaves = cat(1,allSGWaves,sgWaves);
        fgWaves = {dirFG(ovTrials(i),overlappingPts)};
        allFGWaves = cat(1,allFGWaves,fgWaves);
    end
end 
allOvWaves = cat(2,allSGWaves,allFGWaves);

m1_10deg = {uniqueDirs,allUniqueDirs,allOvWaves,ovTrials};    
% get circular correlation for M1
[rho1(3), pval1(3)] = circ_corrcc(allUniqueDirs(1,:),allUniqueDirs(2,:));

clear uniqueDirs uniqueDirsTemp allFGWaves allSGWaves allUniqueDirs dirFG overlappingPts dirSG emptyCells fgWaves i j numTrials outputsM1 ovTrials sgWaves waveTemp waveVector
% for M2
numTrials = numel(outputsM2(1,:)); % same number of trials in both output structures

%for 0 deg
%initialize outputs
waveVector = nan(numTrials,length(timeVals),numFrequencyRanges);
uniqueDirs = cell(numFrequencyRanges,numTrials);
waveBounds = cell(numFrequencyRanges,numTrials);

for i = 1:numTrials
    for j = 1:numFrequencyRanges
    [waveVector(i,:,j),uniqueDirs{j,i},waveBounds{j,i}] = getWaveSegments(outputsM2{j,i},timeVals,wobbleLim1,segOption1,boundryLims, lengthLimit);
    end
end

% find overlapping waves 
allUniqueDirs = [];
dirSG = nan(numTrials,length(timeVals));
dirFG = nan(numTrials,length(timeVals));
waveBoundsOv = cell(1,numTrials);
emptyCells = zeros(1,numTrials);
overlap = 0.5;
for i = 1:numTrials
    [waveBoundsOv{1,i},dirSG(i,:),dirFG(i,:),uniqueDirsTemp,emptyCells(i)] = getOverlappingWaves(waveVector(i,:,1),waveBounds{1,i},waveVector(i,:,2),waveBounds{2,i},overlap);   
    allUniqueDirs = cat(2,allUniqueDirs,uniqueDirsTemp);
end
allUniqueDirs(:,isnan(allUniqueDirs(1,:))) = [];

% get all wave angles
waveBoundsOv(emptyCells==1) = [];
ovTrials = find(emptyCells==0);

allSGWaves = [];
allFGWaves = [];
for i = 1:length(waveBoundsOv)
    waveTemp = waveBoundsOv{i};
    for j = 1:size(waveTemp{1},2)
        overlappingPts = intersect(waveTemp{1}(1,j):waveTemp{1}(2,j),waveTemp{2}(1,j):waveTemp{2}(2,j));
        sgWaves = {dirSG(ovTrials(i),overlappingPts)};
        allSGWaves = cat(1,allSGWaves,sgWaves);
        fgWaves = {dirFG(ovTrials(i),overlappingPts)};
        allFGWaves = cat(1,allFGWaves,fgWaves);
    end
end 
allOvWaves = cat(2,allSGWaves,allFGWaves);

m2_0deg = {uniqueDirs,allUniqueDirs,allOvWaves,ovTrials};    
% get circular correlation for M1
[rho2(1), pval2(1)] = circ_corrcc(allUniqueDirs(1,:),allUniqueDirs(2,:));

% for 5 deg
%initialize outputs
waveVector = nan(numTrials,length(timeVals),numFrequencyRanges);
uniqueDirs = cell(numFrequencyRanges,numTrials);
waveBounds = cell(numFrequencyRanges,numTrials);

for i = 1:numTrials
    for j = 1:numFrequencyRanges
    [waveVector(i,:,j),uniqueDirs{j,i},waveBounds{j,i}] = getWaveSegments(outputsM2{j,i},timeVals,wobbleLim2,segOption2,boundryLims, lengthLimit);
    end
end

% find overlapping waves 
allUniqueDirs = [];
dirSG = nan(numTrials,length(timeVals));
dirFG = nan(numTrials,length(timeVals));
waveBoundsOv = cell(1,numTrials);
emptyCells = zeros(1,numTrials);
overlap = 0.5;
for i = 1:numTrials
    [waveBoundsOv{1,i},dirSG(i,:),dirFG(i,:),uniqueDirsTemp,emptyCells(i)] = getOverlappingWaves(waveVector(i,:,1),waveBounds{1,i},waveVector(i,:,2),waveBounds{2,i},overlap);   
    allUniqueDirs = cat(2,allUniqueDirs,uniqueDirsTemp);
end
allUniqueDirs(:,isnan(allUniqueDirs(1,:))) = [];

% get all wave angles
waveBoundsOv(emptyCells==1) = [];
ovTrials = find(emptyCells==0);

allSGWaves = [];
allFGWaves = [];
for i = 1:length(waveBoundsOv)
    waveTemp = waveBoundsOv{i};
    for j = 1:size(waveTemp{1},2)
        overlappingPts = intersect(waveTemp{1}(1,j):waveTemp{1}(2,j),waveTemp{2}(1,j):waveTemp{2}(2,j));
        sgWaves = {dirSG(ovTrials(i),overlappingPts)};
        allSGWaves = cat(1,allSGWaves,sgWaves);
        fgWaves = {dirFG(ovTrials(i),overlappingPts)};
        allFGWaves = cat(1,allFGWaves,fgWaves);
    end
end 
allOvWaves = cat(2,allSGWaves,allFGWaves);

m2_5deg = {uniqueDirs,allUniqueDirs,allOvWaves,ovTrials};    
% get circular correlation for M1
[rho2(2), pval2(2)] = circ_corrcc(allUniqueDirs(1,:),allUniqueDirs(2,:));

% 10 degrees
%initialize outputs
waveVector = nan(numTrials,length(timeVals),numFrequencyRanges);
uniqueDirs = cell(numFrequencyRanges,numTrials);
waveBounds = cell(numFrequencyRanges,numTrials);

for i = 1:numTrials
    for j = 1:numFrequencyRanges
    [waveVector(i,:,j),uniqueDirs{j,i},waveBounds{j,i}] = getWaveSegments(outputsM2{j,i},timeVals,wobbleLim3,segOption2,boundryLims, lengthLimit);
    end
end

% find overlapping waves 
allUniqueDirs = [];
dirSG = nan(numTrials,length(timeVals));
dirFG = nan(numTrials,length(timeVals));
waveBoundsOv = cell(1,numTrials);
emptyCells = zeros(1,numTrials);
overlap = 0.5;
for i = 1:numTrials
    [waveBoundsOv{1,i},dirSG(i,:),dirFG(i,:),uniqueDirsTemp,emptyCells(i)] = getOverlappingWaves(waveVector(i,:,1),waveBounds{1,i},waveVector(i,:,2),waveBounds{2,i},overlap);   
    allUniqueDirs = cat(2,allUniqueDirs,uniqueDirsTemp);
end
allUniqueDirs(:,isnan(allUniqueDirs(1,:))) = [];

% get all wave angles
waveBoundsOv(emptyCells==1) = [];
ovTrials = find(emptyCells==0);

allSGWaves = [];
allFGWaves = [];
for i = 1:length(waveBoundsOv)
    waveTemp = waveBoundsOv{i};
    for j = 1:size(waveTemp{1},2)
        overlappingPts = intersect(waveTemp{1}(1,j):waveTemp{1}(2,j),waveTemp{2}(1,j):waveTemp{2}(2,j));
        sgWaves = {dirSG(ovTrials(i),overlappingPts)};
        allSGWaves = cat(1,allSGWaves,sgWaves);
        fgWaves = {dirFG(ovTrials(i),overlappingPts)};
        allFGWaves = cat(1,allFGWaves,fgWaves);
    end
end 
allOvWaves = cat(2,allSGWaves,allFGWaves);

m2_10deg = {uniqueDirs,allUniqueDirs,allOvWaves,ovTrials};    
% get circular correlation for M1
[rho2(3), pval2(3)] = circ_corrcc(allUniqueDirs(1,:),allUniqueDirs(2,:));
clear uniqueDirs uniqueDirsTemp allFGWaves allSGWaves allUniqueDirs dirFG overlappingPts dirSG emptyCells fgWaves i j numTrials outputsM1 ovTrials sgWaves waveTemp waveVector
%% plot results
binWidth = 10;
tiledlayout(2,5,'TileSpacing','Compact');
colorVals = cat(1,[52 148 186]./255,[236 112 22]./255);
h1 = nexttile;
sgMean = cell2mat(m1_0deg{1}(1,:));
fgMean = cell2mat(m1_0deg{1}(2,:));

makePolarPlot({fgMean,sgMean},binWidth,h1,colorVals,1)
title('Directions - Unsorted')

h2 = nexttile;
sgOvMean = m1_0deg{2}(1,:);
fgOvMean = m1_0deg{2}(2,:);
makePolarPlot({fgOvMean,sgOvMean},binWidth,h2,colorVals,1)
title('Directions - Overlapping')

% get scatter plot - 0 deg
nexttile
scatter(sgOvMean,fgOvMean,'filled','LineWidth',0.01)
% hold on
% sgOvFull = m1_0deg{3}(:,1);
% fgOvFull = m1_0deg{3}(:,2);
randWaves = [1 2 3 4 6];
% % randWaves = randsample(1:length(sgOvFull),5);
% colors = jet(numel(randWaves));
% for i = 1:numel(randWaves)
%     % scatter(sgOvMean(randWaves(i)),fgOvMean(randWaves(i)),'filled','LineWidth',0.01,'MarkerFaceColor',colors(i,:))
%     s = scatter(wrapToPi(sgOvFull{randWaves(i)}),wrapToPi(fgOvFull{randWaves(i)}),'filled','MarkerFaceColor',colors(i,:),'MarkerEdgeColor','none');
%     s.MarkerFaceAlpha = 0.05;
%     hold on
% end
xlabel('Slow Gamma')
ylabel('Fast Gamma')
title('SG/FG Directions')
xticks([-3.14 0 3.14])
xticklabels({'-180','0','180'})
yticks([-3.14 0 3.14])
yticklabels({'-180','0','180'})
annotation('textbox', [0.48, 0.75, 0.1, 0.1], 'String', ['r:',num2str(round(rho1(1),3)),', p:',num2str(round(pval1(1),3))],'EdgeColor','none')
xlim([-3.14 3.14])
ylim([-3.14 3.14])
axis square

% get scatter plot - 5 deg
nexttile
sgOvMean = m1_5deg{2}(1,:);
fgOvMean = m1_5deg{2}(2,:);
scatter(sgOvMean,fgOvMean,'filled','LineWidth',0.01)
hold on
sgOvFull = m1_5deg{3}(:,1);
fgOvFull = m1_5deg{3}(:,2);
% randWaves = randsample(1:length(sgOvFull),5);
% colors = parula(numel(randWaves));
for i = 1:numel(randWaves)
    scatter(wrapToPi(sgOvMean(randWaves(i))),wrapToPi(fgOvMean(randWaves(i))),'filled','LineWidth',0.01,'MarkerFaceColor',colors(i,:))
    hold on
    s = scatter(wrapToPi(sgOvFull{randWaves(i)}),wrapToPi(fgOvFull{randWaves(i)}),'filled','MarkerFaceColor',colors(i,:),'MarkerEdgeColor','none');
    s.MarkerFaceAlpha = 0.05;
end
xlabel('Slow Gamma')
ylabel('Fast Gamma')
title('SG/FG Directions')
xticks([-3.14 0 3.14])
xticklabels({'-180','0','180'})
yticks([-3.14 0 3.14])
yticklabels({'-180','0','180'})
annotation('textbox', [0.65, 0.75, 0.1, 0.1], 'String', ['r:',num2str(round(rho1(2),3)),', p:',num2str(round(pval1(2),3))],'EdgeColor','none')
xlim([-3.14 3.14])
ylim([-3.14 3.14])
axis square

% get scatter plot - 10 deg
nexttile
sgOvMean = m1_10deg{2}(1,:);
fgOvMean = m1_10deg{2}(2,:);
scatter(sgOvMean,fgOvMean,'filled','LineWidth',0.01)
hold on
sgOvFull = m1_10deg{3}(:,1);
fgOvFull = m1_10deg{3}(:,2);
% randWaves = randsample(1:length(sgOvFull),5);
% colors = parula(numel(randWaves));
for i = 1:numel(randWaves)
    scatter(wrapToPi(sgOvMean(randWaves(i))),wrapToPi(fgOvMean(randWaves(i))),'filled','LineWidth',0.01,'MarkerFaceColor',colors(i,:))
    hold on
    s = scatter(wrapToPi(sgOvFull{randWaves(i)}),wrapToPi(fgOvFull{randWaves(i)}),'filled','MarkerFaceColor',colors(i,:),'MarkerEdgeColor','none');
    s.MarkerFaceAlpha = 0.05;
end
xlabel('Slow Gamma')
ylabel('Fast Gamma')
title('SG/FG Directions')
xticks([-3.14 0 3.14])
xticklabels({'-180','0','180'})
yticks([-3.14 0 3.14])
yticklabels({'-180','0','180'})
annotation('textbox', [0.82, 0.75, 0.1, 0.1], 'String', ['r:',num2str(round(rho1(3),3)),', p:',num2str(round(pval1(3),3))],'EdgeColor','none')
xlim([-3.14 3.14])
ylim([-3.14 3.14])
axis square


h3 = nexttile;
sgMean = cell2mat(m2_0deg{1}(1,:));
fgMean = cell2mat(m2_0deg{1}(2,:));

makePolarPlot({fgMean,sgMean},binWidth,h3,colorVals,1)
title('Directions - Unsorted')

h4 = nexttile;
sgOvMean = m2_0deg{2}(1,:);
fgOvMean = m2_0deg{2}(2,:);
makePolarPlot({fgOvMean,sgOvMean},binWidth,h4,colorVals,1)
title('Directions - Overlapping')

% get scatter plot - 0 deg
nexttile
scatter(sgOvMean,fgOvMean,'filled','LineWidth',0.01)
% hold on
% sgOvFull = m2_0deg{3}(:,1);
% fgOvFull = m2_0deg{3}(:,2);
% randWaves = randsample(1:length(sgOvFull),5);
% colors = parula(numel(randWaves));
% for i = 1:numel(randWaves)
%     % scatter(sgOvMean(randWaves(i)),fgOvMean(randWaves(i)),'filled','LineWidth',0.01,'MarkerFaceColor',colors(i,:))
% 
%     s = scatter(wrapToPi(sgOvFull{randWaves(i)}),wrapToPi(fgOvFull{randWaves(i)}),'filled','MarkerFaceColor',colors(i,:),'MarkerEdgeColor','none');
%     s.MarkerFaceAlpha = 0.05;
%  hold on
% end
xlabel('Slow Gamma')
ylabel('Fast Gamma')
title('SG/FG Directions')
xticks([-3.14 0 3.14])
xticklabels({'-180','0','180'})
yticks([-3.14 0 3.14])
yticklabels({'-180','0','180'})
annotation('textbox', [0.46, 0.13, 0.1, 0.1], 'String', ['r:',num2str(round(rho2(1),3)),', p:',num2str(round(pval2(1),3))],'EdgeColor','none')
xlim([-3.14 3.14])
ylim([-3.14 3.14])
axis square

% get scatter plot - 5 deg
nexttile
sgOvMean = m2_5deg{2}(1,:);
fgOvMean = m2_5deg{2}(2,:);
scatter(sgOvMean,fgOvMean,'filled','LineWidth',0.01)
hold on
sgOvFull = m2_5deg{3}(:,1);
fgOvFull = m2_5deg{3}(:,2);
% randWaves = randsample(1:length(sgOvFull),5);
% colors = parula(numel(randWaves));
for i = 1:numel(randWaves)
    scatter(sgOvMean(randWaves(i)),fgOvMean(randWaves(i)),'filled','LineWidth',0.01,'MarkerFaceColor',colors(i,:))
    hold on
    s = scatter(wrapToPi(sgOvFull{randWaves(i)}),wrapToPi(fgOvFull{randWaves(i)}),'filled','MarkerFaceColor',colors(i,:),'MarkerEdgeColor','none');
    s.MarkerFaceAlpha = 0.05;
end
xlabel('Slow Gamma')
ylabel('Fast Gamma')
title('SG/FG Directions')
xticks([-3.14 0 3.14])
xticklabels({'-180','0','180'})
yticks([-3.14 0 3.14])
yticklabels({'-180','0','180'})
annotation('textbox', [0.63, 0.13, 0.1, 0.1], 'String', ['r:',num2str(round(rho2(2),3)),', p:',num2str(round(pval2(2),3))],'EdgeColor','none')
xlim([-3.14 3.14])
ylim([-3.14 3.14])
axis square

% get scatter plot - 10 deg
nexttile
sgOvMean = m2_10deg{2}(1,:);
fgOvMean = m2_10deg{2}(2,:);
scatter(sgOvMean,fgOvMean,'filled','LineWidth',0.01)
hold on
sgOvFull = m2_10deg{3}(:,1);
fgOvFull = m2_10deg{3}(:,2);
% randWaves = randsample(1:length(sgOvFull),5);
% colors = parula(numel(randWaves));
for i = 1:numel(randWaves)
    scatter(sgOvMean(randWaves(i)),fgOvMean(randWaves(i)),'filled','LineWidth',0.01,'MarkerFaceColor',colors(i,:))
    hold on
    s = scatter(wrapToPi(sgOvFull{randWaves(i)}),wrapToPi(fgOvFull{randWaves(i)}),'filled','MarkerFaceColor',colors(i,:),'MarkerEdgeColor','none');
    s.MarkerFaceAlpha = 0.05;
end
xlabel('Slow Gamma')
ylabel('Fast Gamma')
title('SG/FG Directions')
xticks([-3.14 0 3.14])
xticklabels({'-180','0','180'})
yticks([-3.14 0 3.14])
yticklabels({'-180','0','180'})
annotation('textbox', [0.80, 0.13, 0.1, 0.1], 'String', ['r:',num2str(round(rho2(3),3)),', p:',num2str(round(pval2(3),3))],'EdgeColor','none')
xlim([-3.14 3.14])
ylim([-3.14 3.14])
axis square


annotation('textbox', [0.02, 0.565, 0.2, 0.2], 'String', 'M1','FontSize',24,'FontWeight','bold','EdgeColor','none')
annotation('textbox', [0.02, 0.11, 0.2, 0.2], 'String', 'M2','FontSize',24,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.025, 0.95, 0.02, 0.02], 'string', 'A','FontSize',24,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.22, 0.965, 0.02, 0.02], 'string', 'B','FontSize',24,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.4, 0.965, 0.02, 0.02], 'string', 'C','FontSize',24,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.58, 0.965, 0.02, 0.02], 'string', 'D','FontSize',24,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.74, 0.965, 0.02, 0.02], 'string', 'E','FontSize',24,'FontWeight','bold','EdgeColor','none')




annotation('textbox',[0.025, 0.5, 0.02, 0.02], 'string', 'F','FontSize',24,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.22, 0.5, 0.02, 0.02], 'string', 'G','FontSize',24,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.4, 0.5, 0.02, 0.02], 'string', 'H','FontSize',24,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.58, 0.5, 0.02, 0.02], 'string', 'I','FontSize',24,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.74, 0.5, 0.02, 0.02], 'string', 'J','FontSize',24,'FontWeight','bold','EdgeColor','none')
annotation('rectangle',[0.92, 0.75, 0.02, 0.02], 'FaceColor',colorVals(1,:),'EdgeColor','black')
annotation('rectangle',[0.92, 0.7, 0.02, 0.02], 'FaceColor',colorVals(2,:),'EdgeColor','black')
annotation('textbox',[0.94, 0.73, 0.1, 0.05], 'string', 'Slow Gamma','FontSize',10,'EdgeColor','none')
annotation('textbox',[0.94, 0.68, 0.1, 0.05], 'string', 'Fast Gamma','FontSize',10,'EdgeColor','none')
