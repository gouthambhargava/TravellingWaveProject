function [nanScored] = zscoreNan(data) % zscoring when data has nans
meanData = mean(data(~isnan(data)));
stdData = std(data(~isnan(data)));
nanScored = (data-meanData)/stdData;
end