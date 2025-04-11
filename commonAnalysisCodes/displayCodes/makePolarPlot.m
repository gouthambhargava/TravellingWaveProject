function makePolarPlot(allPhases,binWidth,plotAxis,color)
% plot polar histogram and the mean direction
% Inputs - phaseValues: Give phase values as a cell. Multiple phase vectors
% may be given in a cell, in which case they will be plotted in different
% colors
% binWidth: set at 10deg by default
% plotAxis: specify the plot handle or leave empty
% number of colors and number of phase value cells should be equal

if nargin<4
    colorVals = parula(numel(allPhases));
end

if isempty(plotAxis) % create axis
    plotAxis = getPlotHandles(1,1,[0.1 0.1 0.8 0.8],0.01,0.1,0);
end


for phasei = 1:numel(allPhases)
    phaseValues = cell2mat(allPhases{phasei});
    color = colorVals(phasei,:);
    phaseValues = reshape(phaseValues,[1,numel(phaseValues)]);
    phaseValues(isnan(phaseValues)) = [];

    %bin the data and get counts
    nBins = 360/binWidth;
    allBins = linspace(0,360,nBins);
    wrappedAngles = wrapTo360(rad2deg(phaseValues));
    meanAngle = circ_mean(deg2rad(wrappedAngles)');
    meanAngle = deg2rad(mod(rad2deg(meanAngle),360));
    binnedVals = discretize(wrappedAngles,allBins);
    binnedVals(isnan(binnedVals)) = [];
    uniqueVals = unique(binnedVals);
    counts = zeros(1,numel(uniqueVals));
    for i = 1:numel(uniqueVals)
        counts(uniqueVals(i)) = numel(find(binnedVals==uniqueVals(i)));
    end
    counts = counts/max(counts); %normalize counts

    %plot the counts
    %get bin co-ordinates
    binCord = exp(1i*deg2rad([allBins,allBins(2)]));

    for i = 1:numel(counts)
        xcords = counts(i)*[real(binCord(i)),real(binCord(i+1))];
        ycords = counts(i)*[imag(binCord(i)),imag(binCord(i+1))];
        patch([0 xcords 0],[0 ycords 0],color,'parent',plotAxis)
        hold(plotAxis,'on')
    end

    arrowColor = color+(1-color)*0.2;
    quiver(0,0,cos(meanAngle),sin(meanAngle),'AutoScaleFactor',0.9,'LineWidth',3,'Color',arrowColor,'parent',plotAxis)
    hold(plotAxis,'on')

    line([-1 1],[0 0],'Color',[0.811764705882353,0.811764705882353,0.811764705882353],'parent',plotAxis)
    line([0 0],[-1 1],'Color',[0.811764705882353,0.811764705882353,0.811764705882353],'parent',plotAxis)
    circle(0,0,1,[0,0,0],plotAxis);
    circle(0,0,0.5,[0.811764705882353,0.811764705882353,0.811764705882353],plotAxis);
    axis(plotAxis,[-1,1,-1,1])
end
    axis(plotAxis,'square')
    axis(plotAxis,'off')

end

function circle(x,y,r,colorVal,plotAxis)
    hold on
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    plot(plotAxis,xunit, yunit,'Color',colorVal);
    hold off
end
