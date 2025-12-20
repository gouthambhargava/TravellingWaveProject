function outputs = getTWCircParams(phiMat,burstMat,timeVals,locList,electrodeFraction,electrodeChoice,arrayType,waveDetectionMethod,nPerm)
% Inputs
% phiMat - phases of single trial data in electrodes x time points format
% burstMat - time series indicting burst times in electrodes x time points format
% timeVals - vector of time values
% goodElectrodes - list of good electrodes
% locList - array of electrodes x 2 containing the x and y position of each
% for MEA. In case of EEG, supply the chanlocs structure.
% electrode for microelectrodes or provide chanlocs file for EEG
% electrodeFraction - identify as a cluster only when the number of
% electrodes (out of all electrodes) showing a burst exceeds this fraction.
% electrodeChoice - 'all' or 'selected': Choose either all electrodes or only the electrodes showing a burst for TW calculations
% nPerm - when not set to [], TW statistics are computed nPerm times after shuffling the phases.
% arrayType - microelectrode or EEG - string
% waveDetectionMethod - set to 1 for circlinear regression method. Set to 2
% to apply cluster based circLinear regression method

if ~exist('nPerm','var'); nPerm = [];  end

% Parameters
fs = round(1/(timeVals(2)-timeVals(1))); %sampling frequency
timePoints = length(timeVals);
searchTimePts = dsearchn(timeVals',0.25):dsearchn(timeVals',0.75);
numGoodElectrodes = length(locList);
neighbourLimit = round(numGoodElectrodes*0.25);

if strcmp(electrodeChoice,'selected')~=1
    burstMat = ones(numGoodElectrodes,timePoints);
end


burstMat(isnan(burstMat)) = 0;
burstVec = sum(burstMat)/numGoodElectrodes;
burstVec(burstVec>electrodeFraction) = 1;
burstVec(burstVec~=1) = nan;

if ~isempty(nPerm)
    pgdPerm = zeros(timePoints,nPerm);
end


% numElectrodeCutoff = round(electrodeFraction*numel(goodElectrodes));

if strcmpi(arrayType,'Microelectrode')
    elecDist = 400*10^-6; %distance between adjacent electrodes in the array in m.
else
    elecDist = 2; % set average interelectrode distance and reduce the dimensions of the electrode position data to facilitate circlinear regression (only done for 2d array)
    locList = [locList(1:numGoodElectrodes).X;locList(1:numGoodElectrodes).Y;locList(1:numGoodElectrodes).Z]';
    [coeff, score] = pca(locList);
    locList = coeff(1:2,1:2)*score(:,1:2)'; % selects the highest two coefficients
end

%initialize results
pgd = nan(numGoodElectrodes,timePoints);
direction = nan(numGoodElectrodes,timePoints);
% cluster = zeros(1,timePoints);
sFreq = nan(numGoodElectrodes,timePoints);
coh = zeros(timePoints,1);

%% run circ reg after getting significant electrodes from burst detection
for timei = 1:numel(searchTimePts)
    if burstVec(searchTimePts(timei))==1
        phiVals = phiMat(:,searchTimePts(timei));
        burstLocs = burstMat(:,searchTimePts(timei));
        % cluster(timei,1) = numel(elecs);
        coh(searchTimePts(timei),1) = abs(mean(exp(1i* phiVals(burstLocs>0))));
        [direction(:,searchTimePts(timei)), pgd(:,searchTimePts(timei)),sFreq(:,searchTimePts(timei))] = getWaveMetrics(locList,phiVals,burstLocs,waveDetectionMethod,neighbourLimit);
        % do regression analysis on -the polar and linear coordinates
        if ~isempty(nPerm)
            permVar = zeros(length(phiVals),nPerm);
            for perm = 1:nPerm
                permVar(:,perm) = phiVals(randperm(length(circularCord)),:);
            end
            for perm = 1:nPerm
                [~, pgd(:,timei),~] = getWaveMetrics(locList,phaseMat,burstLocs,waveDetectionMethod,neighbourLimit,numElectrodeCutoff);
            end
        end
    end
end

%% clean up the output variables. Removes values in PGD under 0 and also in
% direction, speed, sfreq and clusters
pgd(pgd<0) = nan;
direction(isnan(pgd)) = nan;
sFreq(isnan(pgd)) = nan;

%%  get additional params
outputs.direction = direction+pi;
outputs.Wavelength = (1./sFreq).*elecDist;  % in m/rad
deltaT = 1/fs; % in s
phi_dot = circ_mean(circ_dist(phiMat(:,2:end),phiMat(:,1:end-1)));
tempFreq = (phi_dot/deltaT)/(2*pi); % convert to Hz
tempFreq = [tempFreq,nan];
tempFreq(tempFreq==0) = nan;
outputs.speed = (tempFreq./sFreq)/1000; % convert to m/s
outputs.tempFreq = wrapTo360(rad2deg(tempFreq));
outputs.speed(outputs.speed==inf) = nan;
outputs.pgd = pgd;
% outputs.clusters = cluster; %the number of electrodes involved in the TW
outputs.coh = coh;
outputs.sFreq = sFreq;
outputs.burstVec = burstVec;
outputs.elecNum = size(locList,1);

if ~isempty(nPerm)
    outputs.pgdPerm = pgdPerm;
    outputs.pgdPermP = prctile(pgdPerm',0.99);
end
end