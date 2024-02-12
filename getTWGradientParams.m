function [outputs] = getTWGradientParams(data,goodElectrodes,timeVals,freqs,req)
%% ref for this method - Rubino, D., Robbins, K. & Hatsopoulos, N. Propagating waves mediate information transfer in the motor cortex. Nat Neurosci 9, 1549–1557 (2006). https://doi.org/10.1038/nn1802
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
    phiVals = unwrap(angle(hilbert(filteredSignal')))';
    gridLayout = flipud(fliplr(reshape(1:81,[9,9]))); 
    timePoints = size(data,2);
    pgd = zeros(1,timePoints);
    speed = zeros(1,timePoints);
    direction = zeros(1,timePoints);
    clusters = zeros(1,timePoints);
    phiGrid = nan(size(gridLayout,1),size(gridLayout,2),timePoints);
    burstLocs = zeros(size(gridLayout,1),size(gridLayout,2),timePoints);
    % reshape phase values to grid and select only those electrodes which show
    % gamma activity
    for i = 1:size(phiVals,1)
        [x,y] = find(gridLayout==goodElectrodes(i));
        phiGrid(x,y,:) = phiVals(i,:);
        burstLocs(x,y,:) = burstTS(i,:);
    end
%% gradient method
    timeDeriv = cat(3,zeros(size(phiGrid,1),size(phiGrid,2)),diff(phiGrid,1,3));
    for j = 1:size(phiGrid,3)
        sigElecs = find(burstLocs(:,:,j));
        if isempty(sigElecs)
            pgd(1,j) = 0;
            speed(1,j) = 0;
            direction(1,j) = 0;
            clusters(1,j) = 0;
        else
            % Calculate Phase Gradient Directionality(PGD),direction (radians) and wave speed
            gridInt = phiGrid(:,:,j);
            gridInt(setdiff(gridLayout,sigElecs)) = nan;
            [gradx,grady] = gradient(gridInt);
            pgd(1,j) = get_PGD(gradx,grady);
            speed(1,j) = get_waveSpeed(gradx,grady,timeDeriv(:,:,j),elecDist,fs);
            theta = reshape(atan2(grady,gradx),[1,numel(gradx)]);
            theta(isnan(theta)) = [];
            direction(1,j) = circ_mean(theta');
            clusters(1,j) = numel(sigElecs);
        end
    end
    outputs.pgd = pgd;
    outputs.speed = speed;
    outputs.direction = direction;
    outputs.clusters = clusters;
    outputs.analysis = 'Gradient';    
end

function PGD = get_PGD(gradx,grady)
    
%Rubino, D., Robbins, K. & Hatsopoulos, N. Propagating waves mediate information transfer in the motor cortex. Nat Neurosci 9, 1549–1557 (2006). https://doi.org/10.1038/nn1802
    
    % remove nans
    gradx(isnan(grady)) = nan;
    grady(isnan(gradx)) = nan;

    % numerator of the PGD
    sumx = nanmean(nanmean(gradx));
    sumy = nanmean(nanmean(grady));
    thetaNum = sqrt(sumx^2 + sumy^2);
    
    % denomarator of the PGD
    grad2 = sqrt(gradx.^2 + grady.^2);
    thetaDenom = nanmean(nanmean(grad2));
    
    % output pgd
    PGD = thetaNum/thetaDenom;
end

function speed = get_waveSpeed(gradx,grady,diff,distanceElecs,Fs)
 % specify elec_dist in metres (for ex - 400*10^-6)
 % grady and gradx are gradients of the phases across the grid along the x and y
 % axes and diff is the gradient along time 
 % %Rubino, D., Robbins, K. & Hatsopoulos, N. Propagating waves mediate information transfer in the motor cortex. Nat Neurosci 9, 1549–1557 (2006). https://doi.org/10.1038/nn1802
   
    % numerator of the speed
    speedNum = abs(mean(diff(~isnan(diff))*Fs));
 
    % denominator of the speed  
    grad2 = sqrt((gradx/distanceElecs).^2 + (grady/distanceElecs).^2);
    speedDenom = mean((mean(grad2(~isnan(grad2)))));
    speed = speedNum/speedDenom; % given in m/s
end