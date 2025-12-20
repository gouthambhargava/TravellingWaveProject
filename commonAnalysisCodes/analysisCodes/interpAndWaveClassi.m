function [cc, pval,waveType,clusterCube,directionCube,phaseCube] = interpAndWaveClassi(phases,direction,locList,waveBoundries,clusterReq,waveMethod)
% this function requires the wave toolbox - https://github.com/mullerlab/wave-matlab/tree/master

% inputs - phases - phases values m x n (electrodes x time)
%          directionGrid - wave direction values m x n (electrodes x time)
%          waveBoundries - indices of waves over the trial
%          sourceFlag - if set to 1, source is calculated. Otherwise the
%          centre of the grid is taken as the source
%          arrayType - 'microelectrode' or 'EEG'
%          clusterReq - determines whether the classification will be done
%          on only the estimated cluster or the interpolated direction

% Outputs - cc - correlation
%           p - pval
%           waveType - 1 (expanding or planar waveVector) or 2 (spiral
%           wave) or 3 (complex/has spiral and planar elements)

if nargin<5
    clusterReq = 0;
end
if nargin<6
    waveMethod = 2;
end

[x, y] = meshgrid(1:9,1:9);

%% get the major cluster and calculate the wavetype based on that cluster
[directionCube,clusterCube] = interpolateMissingDirs(direction,waveBoundries,locList,1,waveMethod);

% make gird from phase and direction Values
phaseCube = nan(size(x,1),size(x,2),length(phases));

for gridi = 1:length(locList)
    phaseCube(locList(gridi,1),locList(gridi,2),:) = phases(gridi,:);
end

if waveMethod==1
    cc = [];
    pval = [];
    waveType = ones(1,size(waveBoundries,2));
else
    cc = zeros(1,size(waveBoundries,2));
    pval = zeros(1,size(waveBoundries,2));
    waveType = zeros(1,size(waveBoundries,2));
    corr = zeros(2,size(waveBoundries,2));
    allP = zeros(2,size(waveBoundries,2));
    % sourcePlanar = zeros(size(waveBoundries,2),2);
    % sourceSpiral = zeros(size(waveBoundries,2),2);
    % allSources = nan(size(waveBoundries,2),2);
    for waveNo = 1:size(waveBoundries,2)
        phaseInt = circMeanNan(phaseCube(:,:,waveBoundries(1,waveNo):waveBoundries(2,waveNo)),3);
        dirInt = circMeanNan(directionCube(:,:,waveBoundries(1,waveNo):waveBoundries(2,waveNo)),3);
        if clusterReq == 1
            clusterInt = mean(clusterCube(:,:,waveBoundries(1,waveNo):waveBoundries(2,waveNo)),3);
            dirInt(isnan(clusterInt)) = nan;
        end

        % detect expanding waves
        % calculate divergance -> max value over grid which gives the origin of the
        % wave -> circular linear correlation between phase and distance from
        % origin
        div = divergence(x,y,sin(dirInt),cos(dirInt));
        [sourceX, sourceY] = find(max(abs(div(:))));
        
        % sourcePlanar(waveNo,:) = [sourceX, sourceY];
        [corr(1,waveNo),allP(1,waveNo)] = phase_correlation_distance(phaseInt, [sourceX, sourceY],1);
        % detect spiral waves
        % calculate curl -> max value over grid gives origin of wave -> circular
        % circular correlation between
        cur = curl(x,y,cos(dirInt),sin(dirInt));
        [corr(2,waveNo),allP(2,waveNo),~] = phase_correlation_rotation(phaseInt,cur,[],1);
    end
    % allSources = cat(3,sourcePlanar, sourceSpiral);
    corr(allP<0.05) = 0;
    for i = 1:size(corr,2)
        [val,valIndex] = max(corr(:,i));
        if val~=0
            cc(i) = corr(valIndex,i);
            pval(i) = allP(valIndex,i);
            waveType(i) = valIndex;
            % sources(i,:) = allSources(i,:,valIndex);
        else
            cc(i) = 0;
            pval(i) = 0;
            waveType(i) = 3;
        end
    end
end
end






