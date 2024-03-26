function [outputs] = getTWCircParams(filteredSignal,burstData,timeVals,goodElectrodes,nPerm)
%% Inputs
%%Inputs
% data - filtered single trial data in the electrodes x time points format
% goodElectrodes - list of good electrodes
% timeVals - vector of time values
% burstData - time series of burst locations
%% 
if nargin<5
    nPerm = [];
end
%% Parameters
    elecDist = 400*10^-6; %distance between adjacent electrodes in the array in m
    fs = 1/(timeVals(2)-timeVals(1)); %sampling frequency
%% filter and extract instant phase of the lfp data
    phiMat = angle(hilbert(filteredSignal'))';
    gridLayout = rot90(reshape(1:81,[9,9]),2); %set the grid layout
    timePoints = length(timeVals);
    
    %initialize results
    pgd = zeros(timePoints,1);
    direction = zeros(timePoints,1);
    cluster = zeros(timePoints,1);
    circVmean = zeros(timePoints,1);
    sFreq = zeros(timePoints,1);
    mag = zeros(timePoints,1);
    
    if ~isempty(nPerm)
        pgdPerm = zeros(timePoints,nPerm);
    end
    % get location for each electrode for segmentation into clusters 
    locList = nan(length(goodElectrodes),2);
    for i = 1:numel(goodElectrodes)
        [locList(i,1),locList(i,2)] = find(gridLayout==goodElectrodes(i));
    end
    locList = cat(2,locList,ones(length(locList),1));
    clear i

    % run circ reg after getting significant electrodes from burst
    % detection
        for timei = 1:size(phiMat,2)
            phiGrid = phiMat(:,timei);
            elecs = find(burstData(:,timei));          
                   %skip to next time point if cluster has <4 electrodes
                   if numel(elecs)>3
                        cluster(timei,1) = numel(elecs);
                        circularCord = phiGrid(elecs);
                        mag(timei,:) = sqrt(sum(phiGrid.^2));
                        circVmean(timei,1) = circ_mean(circularCord);
                        linearCord = locList(elecs,:,:);
                        linearCord(:,3,:) = [];
                    % do regression analysis on the polar and linear coordinates
                        [direction(timei,1),sFreq(timei,1),~,~,pgd(timei,1)] = circRegMod(circularCord,linearCord);
                        if ~isempty(nPerm) 
                            for perm = 1:nPerm 
                                permVar = zeros(length(circularCord),nPerm);
                                permVar(:,perm) = circularCord(randperm(length(circularCord)));
                            end
                            for perm = 1:nPerm
                                [~,~,~,~,pgdPerm(timei,perm)] = circRegMod(permVar(:,perm),linearCord);
                            end
                        end
                   end
%               display(['Done with timepoint:',num2str(timei)])
        end
%%  get additional params
    outputs.direction = direction+pi;
    outputs.Wavelength = (1./sFreq).*elecDist;  % in m/rad
    deltaT = 1/fs; % in s
    phi_dot= abs(diff(circVmean))/deltaT; 
    outputs.tempFreq = phi_dot/(2*pi); % convert to Hz
    outputs.speed = outputs.tempFreq./(sFreq(2:end)/elecDist); % convert to m/s
    outputs.pgd = pgd;
    outputs.clusters = cluster; %the number of electrodes involved in the TW
    outputs.mag = mag;
    if ~isempty(nPerm)
    outputs.pgdPerm = prctile(pgdPerm',0.99);
    end
end