function getPolarPlot(phaseValues,binWidth,color)

phaseValues = reshape(phaseValues,[1,numel(phaseValues)]);
phaseValues(isnan(phaseValues)) = [];
% phaseValues(phaseValues==0) =[];
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
    patch([0 xcords 0],[0 ycords 0],color)
    hold on
end

arrowColor = color+(1-color)*0.2;
quiver(0,0,cos(meanAngle),sin(meanAngle),'AutoScaleFactor',0.9,'LineWidth',3,'Color',arrowColor)
hold on

line([-1 1],[0 0],'Color',[0.811764705882353,0.811764705882353,0.811764705882353])
line([0 0],[-1 1],'Color',[0.811764705882353,0.811764705882353,0.811764705882353])
circle(0,0,1,[0,0,0]);
circle(0,0,0.5,[0.811764705882353,0.811764705882353,0.811764705882353]);
axis([-1,1,-1,1])
text(1.01,0,['0',char(176)])
% text(0,1.05,['90',char(176)])
text(-1.2,0,['180',char(176)])
% text(0,-1.05,['-270',char(176)])
axis square
axis off


    function circle(x,y,r,colorVal)
        hold on
        th = 0:pi/50:2*pi;
        xunit = r * cos(th) + x;
        yunit = r * sin(th) + y;
        plot(xunit, yunit,'Color',colorVal);
        hold off
    end
end