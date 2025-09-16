function [corrCC,corrP] = circCorrPermute(alpha1, alpha2, nPerm,figureFlag)
if nargin<4
    figureFlag = 0;
end
% Pearsons correlation with monte carlo permutation 
% h0 - no difference between original CC and permuted CC (p>0.05), i.e. no significant correlation
% h1 - significant correlation present (p<0.05)
% calculate circular corr coef nPerm times to get a corrected p value.
% alpha1 and alpha2 are the two vectors of circular values. 
    corrCC = circ_corrcc(alpha1, alpha2);
    allAlpha = cat(2,alpha1,alpha2);
    alphaLength = numel(alpha1);
    corrCCPerm = zeros(1,nPerm);
    for i = 1:nPerm
        allAlpha = allAlpha(randperm(length(allAlpha)));
        alpha1T = allAlpha(1:alphaLength);
        alpha2T = allAlpha(alphaLength+1:end);
        corrCCPerm(i) = circ_corrcc(alpha1T, alpha2T);
    end
    corrP = sum(abs(corrCCPerm)>=abs(corrCC))/(nPerm+1); % is the original CC more than the chance CC obtained by permuting the data

if figureFlag==1
    histogram(corrCCPerm,'Normalization','count')
    hold on
    xline(corrCC,'LineWidth',2,'Color','red')
    xlabel('Bins')
    ylabel('Counts')
    legend('Permuted Values','Original Value')
end

end
