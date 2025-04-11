%% load all summary stats
% get sorted outputs
outputsTW1 = outputsTWA24;
outputsTW2 = outputsTWK24;
numTrials1 = numel(outputsTW1(1,:));
numTrials2 = numel(outputsTW2(1,:));

[newBoundsO1,sigCells1,ovDirSG1,ovDirFG1,dirsOL1] = findOverlappingWaves(outputsTW1,timeVals,0.5,2);
[allDirSG1,dirsSG1,boundsSG1] = getWaveSegments(outputsTW1(1,:),timeVals,2);
[allDirFG1,dirsFG1,boundsFG1] = getWaveSegments(outputsTW1(2,:),timeVals,2);
[rho1, pval1] = circ_corrcc(cell2mat(dirsOL1(1,:))',cell2mat(dirsOL1(2,:))');


[newBoundsO2,sigCells2,ovDirSG2,ovDirFG2,dirsOL2] = findOverlappingWaves(outputsTW2,timeVals,0.5,2);
[allDirSG2,dirsSG2,boundsSG2] = getWaveSegments(outputsTW2(1,:),timeVals,2);
[allDirFG2,dirsFG2,boundsFG2] = getWaveSegments(outputsTW2(2,:),timeVals,2);
[rho2, pval2] = circ_corrcc(cell2mat(dirsOL2(1,:))',cell2mat(dirsOL2(2,:))');

% allPoints1 = getRawOverlap(newBoundsO1,ovDirSG1,ovDirFG1);
% allPoints2 = getRawOverlap(newBoundsO2,ovDirSG2,ovDirFG2);
%% plot direction
thetaBins = 10;
tiledlayout(2,3,'TileSpacing','Compact');
colorVals = cat(1,[52 148 186]./255,[236 112 22]./255);
nexttile
getPolarPlotDos(cell2mat(dirsSG1),cell2mat(dirsFG1),thetaBins,colorVals(1,:),colorVals(2,:))
title('Directions - Unsorted')
nexttile
getPolarPlotDos(cell2mat(dirsOL1(1,:)),cell2mat(dirsOL1(2,:)),thetaBins,colorVals(1,:),colorVals(2,:))
title('Directions - Overlapping')
% get scatter plot
x = cell2mat(dirsOL1(1,:));
y = cell2mat(dirsOL1(2,:));
% mdl = fitlm(x,y);
% coefs = mdl.Coefficients.Estimate; % 2x1 [intercept; slope]
nexttile
scatter(x,y,'filled','LineWidth',0.01)
% hold on
% refline(coefs(2),coefs(1)); 
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

nexttile
getPolarPlotDos(cell2mat(dirsSG2),cell2mat(dirsFG2),thetaBins,colorVals(1,:),colorVals(2,:))
title('Directions - Unsorted')
nexttile
getPolarPlotDos(cell2mat(dirsOL2(1,:)),cell2mat(dirsOL2(2,:)),thetaBins,colorVals(1,:),colorVals(2,:))
title('Directions - Overlapping')
% get scatter plot
x = cell2mat(dirsOL2(1,:));
y = cell2mat(dirsOL2(2,:));
% mdl = fitlm(x,y);
% coefs = mdl.Coefficients.Estimate; % 2x1 [intercept; slope]

nexttile
scatter(x,y,'filled','LineWidth',0.01)
% hold on
% refline(coefs(2),coefs(1)); 
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

allDirSG1(isnan(allDirSG1)) = []; 
allDirFG1(isnan(allDirFG1)) = [];
allDirSG2(isnan(allDirSG2)) = []; 
allDirFG2(isnan(allDirFG2)) = [];
axes('Position',[0.03,0.8,0.1,0.1])
getPolarPlotDos(allDirSG1,allDirFG1,thetaBins,colorVals(1,:),colorVals(2,:),0)
axes('Position',[0.03,0.35,0.1,0.1])
getPolarPlotDos(allDirSG2,allDirFG2,thetaBins,colorVals(1,:),colorVals(2,:),0)
axes('Position',[0.32,0.8,0.1,0.1])
getPolarPlotDos(allPoints1(1,:),allPoints1(2,:),thetaBins,colorVals(1,:),colorVals(2,:),0)
axes('Position',[0.32,0.35,0.1,0.1])
getPolarPlotDos(allPoints2(1,:),allPoints2(2,:),thetaBins,colorVals(1,:),colorVals(2,:),0)
%%
annotation('textbox', [0.04, 0.565, 0.2, 0.2], 'String', 'M1','FontSize',24,'FontWeight','bold','EdgeColor','none')
annotation('textbox', [0.04, 0.11, 0.2, 0.2], 'String', 'M2','FontSize',24,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.025, 0.95, 0.02, 0.02], 'string', 'A','FontSize',24,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.3, 0.965, 0.02, 0.02], 'string', 'B','FontSize',24,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.65, 0.965, 0.02, 0.02], 'string', 'C','FontSize',24,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.025, 0.5, 0.02, 0.02], 'string', 'D','FontSize',24,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.3, 0.5, 0.02, 0.02], 'string', 'E','FontSize',24,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.65, 0.5, 0.02, 0.02], 'string', 'F','FontSize',24,'FontWeight','bold','EdgeColor','none')
