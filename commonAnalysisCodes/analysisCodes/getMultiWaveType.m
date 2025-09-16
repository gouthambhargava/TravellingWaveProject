function [cc, pval,waveType,sources] = getMultiWaveType(phases,direction,locList,waveBoundries,clusterFlag,arrayType)
% this function requires the wave toolbox - https://github.com/mullerlab/wave-matlab/tree/master

% inputs - phases - phases values m x n (electrodes x time) 
%          directionGrid - wave direction values m x n (electrodes x time)
%          waveBoundries - indices of waves over the trial
%          clusterFlag - if set to 1, source is calculated otherwise the
%          centre of the grid is taken as the source
%          arrayType - 'microelectrode' or 'EEG'

% Outputs - cc - correlation 
%           p - pval
%           waveType - 1 (expanding or planar waveVector) or 2 (spiral
%           wave) or 3 (complex/has spiral and planar elements)
if nargin<6
    arrayType = 'microelectrode';
end

if strcmp(arrayType,'microelectrode')==1 % currently only defined for microelectrode, EEG part to be added
    [x, y] = meshgrid(1:9,1:9);
end

if min(size(direction))==1
    direction = reshape(direction,[1,size(direction)]); % make a row vector
    direction = repmat(direction,size(phases,1),1);
end

% make gird from phase and direction Values
phaseCube = nan(size(x,1),size(x,2),length(phases));
directionCube = nan(size(x,1),size(x,2),length(direction));

for gridi = 1:length(locList)
    phaseCube(locList(gridi,1),locList(gridi,2),:) = phases(gridi,:);
    directionCube(locList(gridi,1),locList(gridi,2),:) = direction(gridi,:);
end 

cc = zeros(1,size(waveBoundries,2));
pval = zeros(1,size(waveBoundries,2));
waveType = zeros(1,size(waveBoundries,2));
corr = zeros(2,size(waveBoundries,2));
allP = zeros(2,size(waveBoundries,2));
sourcePlanar = zeros(size(waveBoundries,2),2);
sourceSpiral = zeros(size(waveBoundries,2),2);
sources = nan(size(waveBoundries,2),2);
for waveNo = 1:size(waveBoundries,2)
    phaseInt = circMeanNan(phaseCube(:,:,waveBoundries(1,waveNo):waveBoundries(2,waveNo)),3);
    dirInt = circMeanNan(directionCube(:,:,waveBoundries(1,waveNo):waveBoundries(2,waveNo)),3);

    %% detect expanding waves
    % calculate divergance -> max value over grid which gives the origin of the
    % wave -> circular linear correlation between phase and distance from
    % origin 
    if clusterFlag==1
        div = divergence(x,y,cos(dirInt),sin(dirInt));
        [sourceX, sourceY] = find(max(abs(div(:))));
    else
        sourceX = 5;
        sourceY = 5;
    end
    sourcePlanar(waveNo,:) = [sourceX, sourceY];
    [corr(1,waveNo),allP(1,waveNo)] = phase_correlation_distance(phaseInt, [sourceX, sourceY],1);
    %% detect spiral waves
    % calculate curl -> max value over grid gives origin of wave -> circular
    % circular correlation between 
    cur = curl(x,y,cos(dirInt),sin(dirInt));
    [corr(2,waveNo),allP(2,waveNo),sourceSpiral(waveNo,:)] = phase_correlation_rotation(phaseInt,cur,[],1);
end
allSources = cat(3,sourcePlanar, sourceSpiral);
corr(allP<0.05) = 0;
for i = 1:size(corr,2)
    [val,valIndex] = max(corr(:,i));
    if val~=0
        cc(i) = corr(valIndex,i);
        pval(i) = allP(valIndex,i);
        waveType(i) = valIndex;
        sources(i,:) = allSources(i,:,valIndex);
    else 
        cc(i) = 0;
        pval(i) = 0;
        waveType(i) = 3;
    end
end
end






