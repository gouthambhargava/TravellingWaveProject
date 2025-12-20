function [waveDirs] = getClusterWaves(waveBounds,allDirs)
% get a matrix of direction values based on segmented wave indices
if ~iscell(waveBounds)
    waveBounds = {waveBounds};
end

waveDirs = nan(allDirs);
for i = 1:numel(waveBounds)
    for j = 1:size(waveBounds{i},2)
        waveDirs(:,waveBounds{i}(1,j):waveBounds{i}(2,j),i) = allDirs(:,waveBounds{i}(1,j):waveBounds{i}(2,j),i);
    end
end
