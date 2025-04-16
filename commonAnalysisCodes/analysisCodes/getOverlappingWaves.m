function [newBounds,allDirSg,allDirFg,allUniqueDirs] = getOverlappingWaves(dirSg,sgBounds,dirFg,fgBounds,overlap)
%Find the time indices of overlapping waves between slow and fast gamma
%waves. 
% Inputs - the vector of direction values for waves identified in slow and fast gamma are given by dirSg and dirFg 
% and the start of end points of each wave is given by boundriesTimesSg and boundriesTimesFg
% Overlap - this defines the amount of overlap between the waves that is to be considered significant. Set to 50% by default.
%overlapOptions - set to 1 - it calculates the overlap between two waves 
if nargin<5
    overlap = 0.5;
end

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
% allOverCounts{i} = overlapTW'; %(find(overlapTW(:,1)),:);

%% find out which times overlap and remove the ones that dont
% sigCells = find(cellfun(@numel,allOverCounts));
% allOverCounts = allOverCounts(sigCells);
% boundriesTimesSg = boundriesTimesSg(sigCells);
% boundriesTimesFg = boundriesTimesFg(sigCells);
% dirSg = dirSg(sigCells,:);
% dirFg = dirFg(sigCells,:);
% % clear i j k fgBounds sgBounds sgCounts fgCounts singleFg singleTW

%% get new boundries, only the TW's which overlap
newBounds = cell(2,size(overlapTW,1));
for j = 1:size(overlapTW,1)
    newBounds{1}(:,j) = sgBounds(:,overlapTW(j,1));
    newBounds{2}(:,j) = fgBounds(:,overlapTW(j,2));
end

% get unique directions of overlapping waves
allDirSg = nan(size(dirSg));
allDirFg = nan(size(dirFg));

uniqueDirs1 = [];
uniqueDirs2 = [];
for j = 1:size(newBounds{1},2)
    bounds1 = newBounds{1}(:,j);
    bounds2 = newBounds{2}(:,j);
    allDirSg(bounds1(1):bounds1(2)) = dirSg(bounds1(1):bounds1(2));
    allDirFg(bounds2(1):bounds2(2)) = dirFg(bounds2(1):bounds2(2));
    vals1 = circ_mean(unique(dirSg(bounds1(1):bounds1(2)))');
    vals2 = circ_mean(unique(dirFg(bounds2(1):bounds2(2)))');
    vals1(isnan(vals1)) = [];
    vals2(isnan(vals2)) = [];
    uniqueDirs1 = cat(2,uniqueDirs1,vals1);
    uniqueDirs2 = cat(2,uniqueDirs2,vals2);
end
allUniqueDirs = cat(1,uniqueDirs1,uniqueDirs2);
end