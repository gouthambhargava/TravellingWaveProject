function [uniqueDirs,rhoVal,pVal,rhoVal2,pVal2] = simpleOvWrapper(outputs,timeVals,wobble,segOption,lengthLimit,overlap, waveMethod,watsonReq)

if nargin <7
    waveMethod =1;
    watsonReq = 0;
    rhoVal2 = nan;
    pVal2 = nan;
end
boundryLims = [0.25 0.75];

if waveMethod == 1
    numTrials = size(outputs,2);
    numFrequencyRanges = size(outputs,1);
    waveVector = nan(numTrials,length(timeVals),numFrequencyRanges);
    waveBounds = cell(numFrequencyRanges,numTrials);
    for i = 1:numTrials
        for j = 1:numFrequencyRanges
            [waveVector(i,:,j),~,waveBounds{j,i}] = getWaveSegments(outputs{j,i},timeVals,wobble,segOption,boundryLims, lengthLimit);
        end
    end

    % find overlapping waves
    uniqueDirs = cell(1,numTrials);
    emptyCells = nan(1,numTrials);
    dirSG = nan(numTrials,length(timeVals));
    dirFG = nan(numTrials,length(timeVals));
    % olBounds = cell(1,numTrials);
    % overlap = 0.5;
    for i = 1:numTrials
        [~,dirSG(i,:),dirFG(i,:),uniqueDirs{i},emptyCells(i)] = getOverlappingWaves(waveVector(i,:,1),waveBounds{1,i},waveVector(i,:,2),waveBounds{2,i},overlap);
    end
    % get circular correlation for M1
    % olBounds(emptyCells==1) =[];
    uniqueDirs(:,emptyCells==1) = [];
    uniqueDirs = cell2mat(uniqueDirs); % this is the average overlapping direction
    % dirSG(isnan(dirFG)) = nan;
    % dirFG(isnan(dirSG)) = nan; % this is all directions in the set

    %[rhoVal, pVal] = circ_corrcc(uniqueDirs(1,:),uniqueDirs(2,:));
    [rhoVal, pVal] = circCorrPermute(uniqueDirs(1,:),uniqueDirs(2,:), 50,0);
    if watsonReq ==1
        % [pVal2,rhoVal2]=watsons_U2_approx_p(uniqueDirs(1,:),uniqueDirs(2,:));
        [pVal2,rhoVal2] = watsons_U2_perm_test(uniqueDirs(1,:),uniqueDirs(2,:),1000);
    end
else
    numTrials = size(outputs,2);
    numFrequencyRanges = size(outputs,1);
    waveVector = nan(numTrials,length(timeVals),numFrequencyRanges);
    waveBounds = cell(numFrequencyRanges,numTrials);
    for i = 1:numTrials
        for j = 1:numFrequencyRanges
            [waveVector(i,:,j),~,waveBounds{j,i}] = getWaveSegments(outputs{j,i},timeVals,wobble,4,boundryLims, 10);
        end
    end

    % find overlapping waves
    emptyCells = nan(1,numTrials);
    dirSG = nan(numTrials,length(timeVals));
    dirFG = nan(numTrials,length(timeVals));
    olBounds = cell(1,numTrials);
    % overlap = 0.5;
    for i = 1:numTrials
        [olBounds{i},dirSG(i,:),dirFG(i,:),~,emptyCells(i)] = getOverlappingWaves(waveVector(i,:,1),waveBounds{1,i},waveVector(i,:,2),waveBounds{2,i},overlap,2);
    end
    olBounds(:,emptyCells==1) = [];
    numTrials = length(olBounds);
    %% for all dirs in the overlapping regions
    allUniqueDirs = cell(2,numTrials);
    for i = 1:numTrials
            dir1 = outputs{1,i}.direction;
            dir2 = outputs{2,i}.direction;
            overlaps = olBounds{i};
            uqDirsTemp1 = [];
            uqDirsTemp2 = [];
            numWaves = size(overlaps{1},2);
            for j = 1:numWaves
                ovTimes = intersect(overlaps{1,1}(1,j):overlaps{1,1}(2,j),overlaps{2,1}(1,j):overlaps{2,1}(2,j));
                sgInt = dir1(:,ovTimes);
                fgInt = dir2(:,ovTimes);

                % case 0 - only consider overlapping time values and take all
                % directions
                sgInt(isnan(fgInt)) = nan;
                fgInt(isnan(sgInt)) = nan;
                sgInt = reshape(sgInt,[1,numel(sgInt)]);
                sgInt(isnan(sgInt)) = [];
                fgInt = reshape(fgInt,[1,numel(fgInt)]);
                fgInt(isnan(fgInt)) = [];

                % in case only the mean directions are being considered, uncomment the following commands
                % case 1 - mean across time and values across electrodes
                % preserved
                % sgInt = circMeanNan(sgInt);
                % fgInt = circMeanNan(FgInt);

                % case 2 - mean across time and electrodes
                % sgInt = circ_mean(reshape(sgInt,[1,numel(sgInt)]));
                % fgInt = circ_mean(reshape(sgInt,[1,numel(sgInt)]));


                uqDirsTemp1 = cat(2,uqDirsTemp1,sgInt);
                uqDirsTemp2 = cat(2,uqDirsTemp2,fgInt);
            end
            allUniqueDirs{1,i} = uqDirsTemp1;
            allUniqueDirs{2,i} = uqDirsTemp2;
    end
    
    %% get CC and pval
    % [rhoVal, pVal] = circCorrPermute(cell2mat(allUniqueDirs(1,:)),cell2mat(allUniqueDirs(2,:)), 100,0);
    [rhoVal, pVal] = circ_corrcc(cell2mat(allUniqueDirs(1,:)),cell2mat(allUniqueDirs(2,:)));
    if watsonReq ==1
        [pVal2,rhoVal2]=watsons_U2_approx_p(cell2mat(allUniqueDirs(1,:)),cell2mat(allUniqueDirs(2,:)));
    end
    uniqueDirs = allUniqueDirs;
end
