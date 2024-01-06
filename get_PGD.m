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