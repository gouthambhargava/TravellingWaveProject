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