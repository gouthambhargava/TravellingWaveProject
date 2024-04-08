function outputs = getTWCircParams(phiMat,burstMat,timeVals,goodElectrodes,locList,electrodeFraction,electrodeChoice,nPerm)
% Inputs
% phiMat - phases of single trial data in electrodes x time points format
% burstMat - time series indicting burst times in electrodes x time points format
% timeVals - vector of time values
% goodElectrodes - list of good electrodes
% locList - array of size electrodes x 2 containing the x and y position of each electrode
% electrodeFraction - identify as a cluster only when the number of
% electrodes showing a burst exceeds this fraction.
% electrodeChoice - 'all' or 'selected': Choose either all electrodes or only the electrodes showing a burst for TW calculations 
% nPerm - when not set to [], TW statistics are computed nPerm times after shuffling the phases.

if ~exist('nPerm','var');               nPerm = [];                     end

% Parameters
elecDist = 400*10^-6; %distance between adjacent electrodes in the array in m
fs = 1/(timeVals(2)-timeVals(1)); %sampling frequency

timePoints = length(timeVals);
numGoodElectrodes = length(goodElectrodes);

%initialize results
pgd = zeros(timePoints,1);
direction = nan(timePoints,1);
cluster = zeros(timePoints,1);
circVmean = zeros(timePoints,1);
sFreq = zeros(timePoints,1);
coh = zeros(timePoints,1);

burstMat(isnan(burstMat)) = 0;

if ~isempty(nPerm)
    pgdPerm = zeros(timePoints,nPerm);
end

numElectrodeCutoff = round(electrodeFraction*numel(goodElectrodes));

% run circ reg after getting significant electrodes from burst detection
for timei = 1:timePoints
    phiGrid = phiMat(:,timei);
    
    if strcmp(electrodeChoice,'all') % take all electrodes
        elecs = 1:numGoodElectrodes;
    else
        elecs = find(burstMat(:,timei));
    end
    
    if numel(elecs)>=numElectrodeCutoff
        cluster(timei,1) = numel(elecs);
        circularCord = phiGrid(elecs);
        circVmean(timei,1) = circ_mean(circularCord);
        coh(timei,1) = abs(mean(exp(1i*circularCord)));
        linearCord = locList(elecs,:);

        % do regression analysis on the polar and linear coordinates
        [direction(timei,1),sFreq(timei,1),~,~,pgd(timei,1)] = circRegMod(circularCord,linearCord);
        if ~isempty(nPerm)
            for perm = 1:nPerm
                permVar = zeros(length(circularCord),2,nPerm);
                permVar(:,:,perm) = linearCord(randperm(length(linearCord)),:);
            end
            parfor perm = 1:nPerm
                [~,~,~,~,pgdPerm(timei,perm)] = circRegMod(circularCord,permVar(:,:,perm));
            end
        end
%         disp(['done with time ',num2str(timeVals(timei))])
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
outputs.coh = coh;

if ~isempty(nPerm)
    outputs.pgdPerm = pgdPerm;
%     outputs.pgdPerm = prctile(pgdPerm',0.99);
end
end