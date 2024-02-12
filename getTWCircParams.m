function [outputs] = getTWCircParams(data,timeVals,goodElectrodes,freqs,req)
%% Inputs
%%Inputs
% data - single trial data in the electrodes x time points format
% goodElectrodes - list of good electrodes
% timeVals - vector of time values
% freqs - vector of frequency limits (ex. [30 60])
% req - 0 (no filtering required, if MP is being used) or 1 (butterworth filter being used)
%% Parameters
    elecDist = 400*10^-6; %distance between adjacent electrodes in the array in m
    fs = 2000; %sampling frequency

%% filter and extract instant phase of the lfp data
    [burstTS, filteredSignal] = runHilbertBurstLength(data,timeVals,freqs,req);
    burstTS(isnan(burstTS)) = 0;
    phiMat = unwrap(angle(hilbert(filteredSignal')))';
    
    gridLayout = flipud(fliplr(reshape(1:81,[9,9]))); %set the grid layout
    timePoints = size(data,2);
    
    %initialize results
    pgd = zeros(timePoints,1);
    speed = zeros(timePoints,1);
    direction = zeros(timePoints,1);
    cluster = zeros(timePoints,1);
    circVmean = zeros(timePoints,1);
    sFreq = zeros(timePoints,1);

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
            circularCord = [];
            linearCord = [];
            phiGrid = phiMat(:,timei);
            elecs = find(burstTS(:,timei));
                   %skip to next time point if cluster has <4 electrodes
                   if numel(elecs)<4
                       direction(timei,1) = nan;
                       sFreq(timei,1) = nan;
                       Rsq(timei,1) = nan;
                       pgd(timei,1) = 0;
                   else
                    cluster(timei,1) = numel(elecs);
                    circularCord = phiGrid(elecs);
                    circVmean(timei,1) = circ_mean(circularCord);
                    linearCord = locList(elecs,:,:);
                    linearCord(:,3,:) = [];
                    % do regression analysis on the polar and linear coordinates
                    [direction(timei,1),sFreq(timei,1),sl,Rsq,pgd(timei,1)] = circRegMod(circularCord,linearCord);
                   end
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
    outputs.analysis = 'Circreg'; 
end
