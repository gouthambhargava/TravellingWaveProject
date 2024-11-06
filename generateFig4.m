%generate fig 4 - final
% get data
% load data for M1
% dataPath = 'G:\monkeyData\data';
% gridType = 'Microelectrode';
% spatial frequency: 0.5 (1), 1(2), 2 (3), 4 (4), 8 (5), all SFs (6). Note that the same code can be used for the size project also later where stimulus size is changed instead of spatial frequency
% 0 (1), 22.5 (2), 45 (3), 67.5 (4), 90 (5), 112.5 (6), 135 (7), 157.5 (8), all orientations (9)
orientations = [0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5];
load('D:\IISC_work\gitScripts\alpaThresh3.mat')
load('D:\IISC_work\gitScripts\kesariThresh3.mat')
outputsTW1 = {outputsTWA21,outputsTWA22,outputsTWA23,outputsTWA24,outputsTWA25,outputsTWA26,outputsTWA27,outputsTWA28};
outputsTW2 = {outputsTWK21,outputsTWK22,outputsTWK23,outputsTWK24,outputsTWK25,outputsTWK26,outputsTWK27,outputsTWK28};
load('D:\IISC_work\gitScripts\timeVals.mat')
clear outputsTWA21 outputsTWA22 outputsTWA23 outputsTWA24 outputsTWA25 outputsTWA26 outputsTWA27 outputsTWA28 outputsTWA14 outputsTWA24 outputsTWA34 outputsTWA44 outputsTWA54
clear outputsTWK21 outputsTWK22 outputsTWK23 outputsTWK24 outputsTWK25 outputsTWK26 outputsTWK27 outputsTWK28 outputsTWK14 outputsTWK24 outputsTWK34 outputsTWK44 outputsTWK54
% load('all_outputs_75t_new_freqs.mat')
%% get different ori results
for i = 1:numel(outputsTW1)
    [~,~,allDirSg,allDirFg,allUniqueDirs] = findOverlappingWaves(outputsTW1{i},timeVals,0.5,2);
    overlapUqDir1{i} = cell2mat(allUniqueDirs);
    [rhoValue(1,i), pValue(1,i)] = circ_corrcc(overlapUqDir1{i}(1,:),overlapUqDir1{i}(2,:));
    allDirSg(isnan(allDirSg)) = [];
    uqDirSG1{i} = wrapToPi(unique(allDirSg));
    allDirFg(isnan(allDirFg)) = [];
    uqDirFG1{i} = wrapToPi(unique(allDirFg));
end
clear i allDirSg allDirFg allUniqueDirs 

for i = 1:numel(outputsTW2)
    [~,~,allDirSg,allDirFg,allUniqueDirs] = findOverlappingWaves(outputsTW2{i},timeVals,0.5,2);
    overlapUqDir2{i} = cell2mat(allUniqueDirs);
    [rhoValue(2,i), pValue(2,i)] = circ_corrcc(overlapUqDir2{i}(1,:),overlapUqDir2{i}(2,:));
    allDirSg(isnan(allDirSg)) = [];
    uqDirSG2{i} = wrapToPi(unique(allDirSg));
    allDirFg(isnan(allDirFg)) = [];
    uqDirFG2{i} = wrapToPi(unique(allDirFg));
end
clear i allDirSg allDirFg allUniqueDirs 
% plot periodograms and box scatter plots with TW counts
tiledlayout(2,3,'TileSpacing','Compact');
colors = parula(numel(orientations));
colorVals = cat(1,[52 148 186]./255,[236 112 22]./255);

legends = {['0' char(176)], ['22.5' char(176)], ['45' char(176)], ['67.5' char(176)], ['90' char(176)], ['112.5' char(176)], ['135' char(176)], ['157.5' char(176)]};
meanReq = 2;
data1 = uqDirSG1;
data2 = uqDirFG1;
data3 = uqDirSG2;
data4 = uqDirFG2;

% plot directions
% violin plots monkey1 slow gamma - ori
nexttile(1)
getViolinPlotPlain(data1,colors,legends,meanReq)
title('Slow Gamma Wave Directions')
ylabel('Directions')
% ylim([-3.14 3.14])
yticks([-3.14,0,3,14])
yticklabels({'-180','0','180'})

% violin plots monkey1 slow gamma - ori
nexttile(2)
getViolinPlotPlain(data2,colors,legends,meanReq)
title('Fast Gamma Wave Directions')
ylabel('Directions')
% ylim([-3.14 3.14])
yticks([-3.14,0,3,14])
yticklabels({'-180','0','180'})

% violin plots monkey2 slow gamma - ori
nexttile(4)
getViolinPlotPlain(data3,colors,legends,meanReq)
title('Slow Gamma Wave Directions')
ylabel('Directions')
% ylim([-3.14 3.14])
yticks([-3.14,0,3,14])
yticklabels({'-180','0','180'})
xlabel('Orientations')

% violin plots monkey2 slow gamma - ori
nexttile(5)
getViolinPlotPlain(data4,colors,legends,meanReq)
title('Fast Gamma Wave Directions')
ylabel('Directions')
% ylim([-3.14 3.14])
yticks([-3.14,0,3,14])
yticklabels({'-180','0','180'})
xlabel('Orientations')


%plot the corrcoef
nexttile(3)
sigPs = find(pValue(1,:)<0.05);
plot(1:numel(orientations),rhoValue(1,:),'-o','LineWidth',1.5,'MarkerIndices',1:numel(orientations),'Color',colorVals(1,:))
if ~isempty(sigPs)
    hold on
    sigPlot = ones(1,length(sigPs));
    sigPlot(sigPs) = 1;
%     sigPlot(sigPs1) = rhoValue(1,sigPs1)+rhoValue(1,sigPs1).*0.5;
    plot(1:numel(sigPlot),sigPlot*0.3,'*','Color','k','LineWidth',1.5)
end
ylim([-1 1])
xticks([1 2 3 4 5 6 7 8])
xticklabels(legends)
xlabel('Orientations')
ylabel('Circular CC')
title('Circ Corr across orientations:M1')

nexttile(6)
sigPs = find(pValue(2,:)<0.05);
plot(1:numel(orientations),rhoValue(2,:),'-o','LineWidth',1.5,'MarkerIndices',1:numel(orientations),'Color',colorVals(1,:))
if ~isempty(sigPs)
    hold on
    sigPlot = ones(1,length(sigPs));
    sigPlot(sigPs) = 1;
%     sigPlot(sigPs1) = rhoValue(1,sigPs1)+rhoValue(1,sigPs1).*0.5;
    plot(1:numel(sigPlot),sigPlot*0.3,'*','Color','k','LineWidth',1.5)
end
ylim([-1 1])
xticks([1 2 3 4 5 6 7 8])
xticklabels(legends)
xlabel('Orientations')
ylabel('Circular CC')
title('Circ Corr across orientations:M2')
