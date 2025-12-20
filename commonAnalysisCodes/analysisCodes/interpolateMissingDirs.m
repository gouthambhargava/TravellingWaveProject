function [fullDirCube,fullClusterCube] = interpolateMissingDirs(dirs,boundries,locList,methodFlag, waveMethod)
% function interpolates wave directions for missing electrodes.
% methodFlag - 1: interpolation using all electrodes arranged as a vector
%              2: interpolation using a grid
% boundries - list of wave indices
%%
if nargin<5; waveMethod=2; end
if waveMethod==2

    gridLayout = rot90(reshape(1:81,[9,9]),2); %set the grid layout
    fullDirCube = nan(size(gridLayout,1),size(gridLayout,2),size(dirs,2));
    fullClusterCube = nan(size(gridLayout,1),size(gridLayout,2),size(dirs,2));

    for waveId = 1:size(boundries,2)
        directions = dirs(:,boundries(1,waveId):boundries(2,waveId));

        % remove sudden jumps in dirs
        vals = ischange(wrapToPi(directions),2);
        directions(vals) = nan;

        % remove dir values in successive time points that are less than 3 time
        % points
        directionsNew = nan(size(directions,1),size(directions,2)+2);
        for i = 1:size(directions,1)
            directionsNew(i,:) = getDataSegments([nan,directions(i,:),nan],3);
        end
        directions = directionsNew(:,2:end-1);

        directionCube = nan(size(gridLayout,1),size(gridLayout,2),size(directions,2));

        for gridi = 1:size(locList,1)
            directionCube(locList(gridi,1),locList(gridi,2),:) = directions(gridi,:);
        end

        directionCubeMean = circMeanNan(directionCube, 3);
        directionCubeMean(directionCubeMean~=0) = 1;

        CC = bwconncomp(directionCubeMean);
        largestCluster = CC.PixelIdxList{1,1};
        clusterId = CC.PixelIdxList;
        tinyClusters = cell2mat(clusterId(cellfun(@numel,clusterId)<2));
        gridLocations = [];
        if ~isempty(tinyClusters)
            for i = 1:numel(tinyClusters)
                [gridLocations(i,1),gridLocations(i,2)] = find(gridLayout==tinyClusters(i));
            end
        end
        directionCube(gridLocations,:) = nan;

        clusterCube = nan(size(gridLayout));
        clusterCube(largestCluster) = 1;

        if methodFlag==1
            % interpolate for missing direction values
            newDirCube = nan(size(directionCube));
            for i = 1:size(directionCube,3)
                dirFrame = directionCube(:,:,i);
                dirFrame = reshape(dirFrame,[numel(dirFrame),1]);
                dirFrame = unwrap(rad2deg(wrapToPi(dirFrame))*pi/180)*180/pi;
                newFrame = fillmissing(dirFrame,'linear');
                newFrame=mod(newFrame,360);
                newFrame(newFrame>180)=newFrame(newFrame>180)-360;
                newFrame = deg2rad(newFrame);
                newFrame = reshape(newFrame,[9,9]);
                newDirCube(:,:,i) = newFrame;
            end
        else
            % interpolate for missing direction values
            newDirCube = nan(size(directionCube));
            for i = 1:size(directionCube,3)
                dirFrame = directionCube(:,:,i);
                dirFrame = unwrap(rad2deg(wrapToPi(dirFrame))*pi/180)*180/pi;
                newFrame = fillmissing(dirFrame,'knn',75,2);
                newFrame=mod(newFrame,360);
                newFrame(newFrame>180)=newFrame(newFrame>180)-360;
                newFrame = deg2rad(newFrame);
                newDirCube(:,:,i) = newFrame;
            end
        end
        fullDirCube(:,:,boundries(1,waveId):boundries(2,waveId)) = newDirCube;
        fullClusterCube(:,:,boundries(1,waveId):boundries(2,waveId)) = repmat(clusterCube,1,1,numel(boundries(1,waveId):boundries(2,waveId)));
    end

else
    gridLayout = rot90(reshape(1:81,[9,9]),2); %set the grid layout
    fullDirCube = nan(size(gridLayout,1),size(gridLayout,2),size(dirs,2));
    fullClusterCube = nan(size(gridLayout,1),size(gridLayout,2),size(dirs,2));

    for gridi = 1:size(locList,1)
        fullDirCube(locList(gridi,1),locList(gridi,2),:) = dirs(gridi,:);
        fullClusterCube(locList(gridi,1),locList(gridi,2),:) = ones(size(dirs(gridi,:)));
    end
end

end

function [waveVector,boundries] = getDataSegments(dataOrig,lengthLimit)

dataOrig = reshape(dataOrig,[length(dataOrig),1]);%convert to column vector
waveVector = nan(1,length(dataOrig));
data = dataOrig;
data(isnan(data)) = 0;
data(data~=0) = 1;
data = find(data==0);
dataEpochs = find(diff(data)>1);
boundries = [data(dataEpochs)+1,data(dataEpochs+1)-1]';
boundries(:,diff(boundries)<lengthLimit) = [];
for k = 1:size(boundries,2)
    waveVector(boundries(1,k):boundries(2,k)) = dataOrig(boundries(1,k):boundries(2,k));
end
end