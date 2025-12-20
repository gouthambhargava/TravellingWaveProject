function [waveVector,uniqueDirs,waveBounds,clusterCube,waveType,directionCube,phaseCube] = getWaveParameters(outputs,phases,locList,timeVals,lengthLimit, wobble, segMethod,boundryLims,waveMethod)
% a script to nest the wave segmentation and wave type functions for easier display codes
% script calculates the wave segments, interpolation (only for cluster
% based method) and gives the direction, phase and clusters for plotting
% gridLayout = rot90(reshape(1:81,[9,9]),2); %set the grid layout
direction = outputs.direction;
waveVector = nan(3,length(timeVals));
if waveMethod == 2
    [~,~,waveBounds] = getWaveSegments(outputs,timeVals,wobble,segMethod,boundryLims,lengthLimit);
    if  ~isempty(waveBounds)
        [~, ~,waveType,clusterCube,directionCube,phaseCube] = interpAndWaveClassi(phases,direction,locList,waveBounds,1,waveMethod);
        waveVector = nan(numel(waveType),length(timeVals));
        % uniqueDirsExp = cell(1,3);
        uniqueDirs = cell(1,numel(waveType));
        for j = 1:size(waveBounds,2)
            waveVector(j,waveBounds(1,j):waveBounds(2,j)) = 1;
            tempUq = circMeanNan(direction(:,waveBounds(1,j):waveBounds(2,j)),2);
            uniqueDirs{1,j} = tempUq;
        end
    else
        clusterCube = [];
    end
else
    [~,uniqueDirs,waveBounds] = getWaveSegments(outputs,timeVals,wobble,segMethod, boundryLims,lengthLimit);
    uniqueDirs = num2cell(uniqueDirs);
    if  ~isempty(waveBounds)
        [~, ~,waveType,clusterCube,directionCube,phaseCube] = interpAndWaveClassi(phases,direction,locList,waveBounds,0,waveMethod);
    end
    waveVector = nan(numel(waveType),length(timeVals));
    for j = 1:size(waveBounds,2)
        waveVector(j,waveBounds(1,j):waveBounds(2,j)) = 1;
    end
end