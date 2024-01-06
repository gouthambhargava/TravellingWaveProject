function [outputs] = getTWParams(data,goodElectrodes,freqs,req,electrodeList)
if nargin<3
    req = 1;
    electrodeList =[];
end
%% Inputs
% load data and filter from loadLFPData
% permReq - 0 (default)/1, permutation testing will be omitted if permReq is set to 0.
% req - 1/2 follows gradient or circular linear regression respectively
% electrodeList - list the subset of electrodes that are to be taken up for
% further analysis, only applicable for linear reg and cir-linear reg
% method
%freqs - frequency limits of filter used. Ex [30 80];
%% assume parameters
elecDist = 400*10^-6; %distance between adjacent electrodes in the array in m
fs = 2000; %sampling frequency
%% filter and extract instant phase of the lfp data
gridLayout = flipud(fliplr(reshape(1:81,[9,9]))); %set the grid layout for alpaH
numTrials = size(data,3);
timePoints = size(data,2);
pgd = zeros(timePoints,numTrials);
speed = zeros(timePoints,numTrials);
direction = zeros(timePoints,numTrials);
phiGrid = nan(size(gridLayout,1),size(gridLayout,2),timePoints);

% get location for each electrode for segmentation in clusters for method 2
% and 3 according to inst freq
locList = nan(numel(gridLayout),2);
for i = 1:numel(goodElectrodes)
    [locList(goodElectrodes(i),1),locList(goodElectrodes(i),2)] = find(gridLayout==goodElectrodes(i));
end
locList = cat(2,locList,ones(length(locList),1));
clear i


%% gradient method
switch req
    case 1 
   %ref for this method - Rubino, D., Robbins, K. & Hatsopoulos, N. Propagating waves mediate information transfer in the motor cortex. Nat Neurosci 9, 1549–1557 (2006). https://doi.org/10.1038/nn1802
tic
for n = 1:length(numTrials)
    for i = 1:numel(goodElectrodes)
        [x,y] = find(gridLayout==goodElectrodes(i));
        phiGrid(x,y,:) = unwrap(angle(hilbert(data(i,:,n))));
    end 
    timeDeriv = cat(3,zeros(size(phiGrid,1),size(phiGrid,2)),diff(phiGrid,1,3));
    for j = 1:size(phiGrid,3)
        % Calculate Phase Gradient Directionality(PGD) and other parameters
        [gradx,grady] = gradient(phiGrid(:,:,j));
        pgd(j,n) = get_PGD(gradx,grady);
        speed(j,n) = get_waveSpeed(gradx,grady,timeDeriv(:,:,j),elecDist,fs);
        theta = reshape(atan2(grady,gradx),[1,numel(gradx)]);
        theta(isnan(theta)) = [];
        direction(j,n) = circ_mean(theta');
    end
end
    toc
    outputs.pgd = pgd;
    outputs.speed = speed;
    outputs.direction = direction;
    outputs.analysis = 'Gradient';    
    %% circular linear reg
    case 2
%% Loop through the trails and times
tic
for triali = 1:size(data,3)
    instFreqMat = nan(size(gridLayout,1),size(gridLayout,2),timePoints-1);
    instphiMat = nan(size(gridLayout,1),size(gridLayout,2),timePoints);
        for i = 1:numel(goodElectrodes)
            [x,y] = find(gridLayout==goodElectrodes(i));
            instphiMat(x,y,:) = unwrap(angle(hilbert(data(i,:,triali))));
            instFreqMat(x,y,:) = instfreq(data(i,:,triali)',fs,'Method','hilbert');
        end 
    % get electrodes for each time point using instFreq
        for timei = 1:size(instFreqMat,3)
            circularCord = [];
            linearCord = [];
            phiGrid = instphiMat(:,:,timei);
            elecs = find(instFreqMat(:,:,timei)>freqs(1) & instFreqMat(:,:,timei)>freqs(2));
% method for defining a cluster of electrodes that are contiguous, set by
% estimating a a group of electrodes seperated by a minimum distance
% (mindist) from each other. Electrodes away from the main cluster will be
% ignored. Analysis is considered for that time point only if there are atleast 4 electrodes in
% a cluster. Present pipeline considers all the electrodes in given by inst
% freq.
%             distMat = squareform(pdist(locList(gridLayout(elecs),:)));
%             distMat(isnan(distMat)) = 0;
%             distMat(distMat>minDist) = 0; % removes electrodes that do
%             not belong in a continuous cluster. To be refined.
%             distMat(distMat>0) = 1;
%             [elecClus clusSize] = conncomp(graph(distMat,'upper'));
%             numClus = find(clusSize==max(clusSize));
%                 if max(clusSize)>4
%                     if numel(numClus)>1
%                         numClus = numClus(1);
%                     end    
%                     clusters(triali,timei) = numel(gridLayout(elecs(find(elecClus==(numClus)))));
%                     linearCord = locList(gridLayout(elecs(find(elecClus==(numClus)))),:,:);
%                     linearCord(:,3,:) = [];
%                     phiMat = instPhiMat(:,:,timei);
%                     circularCord = phiMat(elecs(find(elecClus==(numClus))));
                  
                   %skip to next time point if cluster has <4 electrodes
                   if numel(elecs)<4
                       direction(triali,timei) = nan;
                       sFreq(triali,timei) = nan;
                       Rsq(triali,timei) = nan;
                       pgd(triali,timei) = nan;
                   else
                    cluster(triali,timei) = numel(elecs);
                    circularCord = phiGrid(elecs);
                    circVmean(triali,timei) = circ_mean(circularCord);
                    linearCord = locList(gridLayout(elecs),:,:);
                    linearCord(:,3,:) = [];
                    % do circlinear reg analysis on the coordinates
                    [direction(triali,timei),sFreq(triali,timei),sl,Rsq,pgd(timei,triali)] = circRegMod(circularCord,linearCord);
                   end
        end
end
toc
%  get additional params
    outputs.direction = direction+pi;
    outputs.Wavelength = (1./sFreq).*elecDist;  % in m/rad
    deltaT = 1/fs; % in s
    phi_dot= abs(diff(circVmean,1,2))/deltaT; 
    outputs.tempFreq = phi_dot/(2*pi); % convert to Hz
    outputs.speed = outputs.tempFreq./(sFreq(:,2:end)/elecDist); % convert to m/s
    outputs.pgd = pgd;
    outputs.clusters = cluster; %the number of electrodes involved in the TW 
    outputs.analysis = 'Circreg'; 
end
end
