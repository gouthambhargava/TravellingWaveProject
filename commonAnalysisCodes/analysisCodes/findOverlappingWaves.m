function [newBounds,sigCells,allDirSg,allDirFg,allUniqueDirs] = findOverlappingWaves(outputsTW,timeVals,overlap,wobbleReq)
if nargin<4
    wobbleReq=2;
end
% get overlapping TW's in all trials
numTrials = numel(outputsTW(1,:));

if wobbleReq ==1
    %for segmentation without wobble
[dirSg,~,boundriesTimesSg] = getWaveSegments(outputsTW(1,:),timeVals,1);
[dirFg,~,boundriesTimesFg] = getWaveSegments(outputsTW(2,:),timeVals,1);
elseif wobbleReq ==2
    % for the original segmentation method 
[dirSg,~,boundriesTimesSg] = getWaveSegments(outputsTW(1,:),timeVals,2);
[dirFg,~,boundriesTimesFg] = getWaveSegments(outputsTW(2,:),timeVals,2);
else
    %for newest 5 degree criteria
[dirSg,~,boundriesTimesSg] = getWaveSegments(outputsTW(1,:),timeVals,3);
[dirFg,~,boundriesTimesFg] = getWaveSegments(outputsTW(2,:),timeVals,3);
end

selectedTrials = [find(cellfun(@numel,boundriesTimesSg)==0),find(cellfun(@numel,boundriesTimesFg)==0)];
selectedTrials = unique(selectedTrials);

% remove trials with no bursts
boundriesTimesSg(selectedTrials) = [];
boundriesTimesFg(selectedTrials) = [];
dirSg(selectedTrials,:) = [];
dirFg(selectedTrials,:) = [];
selectedTrials = setdiff(1:numTrials,selectedTrials);
for i = 1:numel(boundriesTimesSg)
    sgBounds = boundriesTimesSg{i};
    fgBounds = boundriesTimesFg{i};
%     sgCounts = zeros(1,size(sgBounds,2));
    fgCounts = zeros(1,size(fgBounds,2));
   
    overlapTW = [];
    for j = 1:size(sgBounds,2)
        singleSg = sgBounds(1,j):sgBounds(2,j);
        sgLength = length(singleSg);
        for k = 1:length(fgCounts)
            singleFg = fgBounds(1,k):fgBounds(2,k);
            fgLength = length(singleFg);
            intersectLength = numel(intersect(singleSg,singleFg));
            if intersectLength>overlap*sgLength && intersectLength>overlap*fgLength
               overlapTW = cat(1,overlapTW,[j,k]);
            end          
        end
    end
    allOverCounts{i} = overlapTW'; %(find(overlapTW(:,1)),:);
end
%% find out which times overlap and remove the ones that dont
sigCells = find(cellfun(@numel,allOverCounts));
allOverCounts = allOverCounts(sigCells);
boundriesTimesSg = boundriesTimesSg(sigCells);
boundriesTimesFg = boundriesTimesFg(sigCells);
dirSg = dirSg(sigCells,:);
dirFg = dirFg(sigCells,:);
% clear i j k fgBounds sgBounds sgCounts fgCounts singleFg singleTW

%% get new boundries, only the TW's which overlap
% for slow gamma
newBounds = cell(2,numel(allOverCounts));
for i = 1:numel(boundriesTimesSg)
    for j = 1:size(allOverCounts{i},2)
        newBounds{1,i}(:,j) = boundriesTimesSg{i}(:,allOverCounts{i}(1,j));
        newBounds{2,i}(:,j) = boundriesTimesFg{i}(:,allOverCounts{i}(2,j));
    end
end

allDirSg = nan(size(dirSg));
allDirFg = nan(size(dirFg));

for i = 1:numel(newBounds(1,:))
    uniqueDirs1 = [];
    uniqueDirs2 = [];
    for j = 1:size(newBounds{1,i},2)
        bounds1 = newBounds{1,i}(:,j);
        bounds2 = newBounds{2,i}(:,j);
        allDirSg(i,bounds1(1):bounds1(2)) = dirSg(i,bounds1(1):bounds1(2));
        allDirFg(i,bounds2(1):bounds2(2)) = dirFg(i,bounds2(1):bounds2(2));
        vals1 = circ_mean(unique(dirSg(i,bounds1(1):bounds1(2)))');
        vals2 = circ_mean(unique(dirFg(i,bounds2(1):bounds2(2)))');
        vals1(isnan(vals1)) = [];
        vals2(isnan(vals2)) = [];
        uniqueDirs1 = cat(2,uniqueDirs1,vals1);
        uniqueDirs2 = cat(2,uniqueDirs2,vals2);
    end
    allUniqueDirs{1,i} = uniqueDirs1;
    allUniqueDirs{2,i} = uniqueDirs2;
end
end