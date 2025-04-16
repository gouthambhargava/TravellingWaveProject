function outputs = getTWCircParams(phiMat,burstMat,timeVals,goodElectrodes,locList,electrodeFraction,electrodeChoice,nPerm)
% Inputs
% phiMat - phases of single trial data in electrodes x time points format
% burstMat - time series indicting burst times in electrodes x time points format
% timeVals - vector of time values
% goodElectrodes - list of good electrodes
% locList - array of size electrodes x 2 containing the x and y position of each electrode
% electrodeFraction - identify as a cluster only when the number of
% electrodes (out of all electrodes) showing a burst exceeds this fraction.
% electrodeChoice - 'all' or 'selected': Choose either all electrodes or only the electrodes showing a burst for TW calculations 
% nPerm - when not set to [], TW statistics are computed nPerm times after shuffling the phases.

if ~exist('nPerm','var'); nPerm = [];  end

% Parameters
elecDist = 400*10^-6; %distance between adjacent electrodes in the array in m
fs = round(1/(timeVals(2)-timeVals(1))); %sampling frequency
timePoints = length(timeVals);
numGoodElectrodes = length(goodElectrodes);

%initialize results
pgd = zeros(timePoints,1);
direction = nan(timePoints,1);
cluster = zeros(timePoints,1);
sFreq = zeros(timePoints,1);
coh = zeros(timePoints,1);

burstMat(isnan(burstMat)) = 0;
burstVec = sum(burstMat);
burstVec(burstVec>0) = 1;
burstVec(burstVec==0) = nan;

if ~isempty(nPerm)
    pgdPerm = zeros(timePoints,nPerm);
end

numElectrodeCutoff = round(electrodeFraction*numel(goodElectrodes));

% run circ reg after getting significant electrodes from burst detection
for timei = 1:timePoints
    phiGrid = phiMat(:,timei);
    
    if strcmpi(electrodeChoice,'all') % take all electrodes
        elecs = 1:numGoodElectrodes;
    else
        elecs = find(burstMat(:,timei));
    end
    
    if numel(elecs)>=numElectrodeCutoff
        cluster(timei,1) = numel(elecs);
        circularCord = phiGrid(elecs);
        coh(timei,1) = abs(mean(exp(1i*circularCord)));
        linearCord = locList(elecs,:);
        if numel(elecs)>3
            % do regression analysis on the polar and linear coordinates
            [direction(timei,1),sFreq(timei,1),~,~,pgd(timei,1)] = circRegMod(circularCord,linearCord);
            if ~isempty(nPerm)
                permVar = zeros(length(circularCord),nPerm);
                for perm = 1:nPerm
                    permVar(:,perm) = circularCord(randperm(length(circularCord)),:);
                end
                parfor perm = 1:nPerm
                    [~,~,~,~,pgdPerm(timei,perm)] = circRegMod(permVar(:,perm),linearCord);
                end
            end
        end
%          disp(['done with time ',num2str(timeVals(timei))])
    end
end

% clean up the output variables. Removes values in PGD under 0 and also in
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
outputs.speed = (tempFreq./sFreq')/1000; % convert to m/s
outputs.tempFreq = wrapTo360(rad2deg(tempFreq));
outputs.speed(outputs.speed==inf) = nan;
outputs.pgd = pgd;
outputs.clusters = cluster; %the number of electrodes involved in the TW
outputs.coh = coh;
outputs.sFreq = sFreq;
outputs.burstVec = burstVec;

if ~isempty(nPerm)
     outputs.pgdPerm = pgdPerm;
     outputs.pgdPermP = prctile(pgdPerm',0.99);
end
end