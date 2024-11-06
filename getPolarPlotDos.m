function getPolarPlotDos(phaseValues1,phaseValues2,binWidth,color1,color2,axisOptions)
if nargin<6
    axisOptions = 1;
end

%for phase values 1
phaseValues1 = reshape(phaseValues1,[1,numel(phaseValues1)]);
phaseValues1(isnan(phaseValues1)) = [];
%bin the data and get counts
nBins = 360/binWidth;
allBins = linspace(0,360,nBins);
wrappedAngles1 = wrapTo360(rad2deg(phaseValues1));
meanAngle1 = circ_mean(deg2rad(wrappedAngles1)');
meanAngle1 = deg2rad(mod(rad2deg(meanAngle1),360));
binnedVals1 = discretize(wrappedAngles1,allBins);
binnedVals1(isnan(binnedVals1)) = [];
uniqueVals1 = unique(binnedVals1);

counts = zeros(1,numel(uniqueVals1));
for i = 1:numel(uniqueVals1)
    counts(uniqueVals1(i)) = numel(find(binnedVals1==uniqueVals1(i)));
end
counts = counts/max(counts); %normalize counts 

%plot the counts
%get bin co-ordinates
binCord = exp(1i*deg2rad([allBins,allBins(2)]));

for i = 1:numel(counts)
    xcords = counts(i)*[real(binCord(i)),real(binCord(i+1))];
    ycords = counts(i)*[imag(binCord(i)),imag(binCord(i+1))];
    patch([0 xcords 0],[0 ycords 0],color1)
    hold on
end



%for phase values 2
phaseValues2 = reshape(phaseValues2,[1,numel(phaseValues2)]);
phaseValues2(isnan(phaseValues2)) = [];
%bin the data and get counts
nBins = 360/binWidth;
allBins = linspace(0,360,nBins);
wrappedAngles2 = wrapTo360(rad2deg(phaseValues2));
meanAngle2 = circ_mean(deg2rad(wrappedAngles2)');
meanAngle2 = deg2rad(mod(rad2deg(meanAngle2),360));
binnedVals2 = discretize(wrappedAngles2,allBins);
binnedVals2(isnan(binnedVals2)) = [];
uniqueVals2 = unique(binnedVals2);


counts = zeros(1,numel(uniqueVals2));
for i = 1:numel(uniqueVals2)
    counts(uniqueVals2(i)) = numel(find(binnedVals2==uniqueVals2(i)));
end
counts = counts/max(counts); %normalize counts 

%plot the counts
%get bin co-ordinates
binCord = exp(1i*deg2rad([allBins,allBins(2)]));

for i = 1:numel(counts)
    xcords = counts(i)*[real(binCord(i)),real(binCord(i+1))];
    ycords = counts(i)*[imag(binCord(i)),imag(binCord(i+1))];
    patch([0 xcords 0],[0 ycords 0],color2)
    hold on
end


arrowColor1 = color1+(1-color1)*0.2;
quiver(0,0,cos(meanAngle1),sin(meanAngle1),'AutoScaleFactor',0.9,'LineWidth',3,'Color',arrowColor1)
hold on
arrowColor2 = color2+(1-color2)*0.2;
quiver(0,0,cos(meanAngle2),sin(meanAngle2),'AutoScaleFactor',0.9,'LineWidth',3,'Color',arrowColor2)
hold on




line([-1 1],[0 0],'Color',[0.811764705882353,0.811764705882353,0.811764705882353])
line([0 0],[-1 1],'Color',[0.811764705882353,0.811764705882353,0.811764705882353])
circle(0,0,1,[0,0,0]);
circle(0,0,0.5,[0.811764705882353,0.811764705882353,0.811764705882353]);
axis([-1,1,-1,1])
if axisOptions == 1
    text(1.01,0,['0',char(176)])
    % text(0,1.05,['90',char(176)])
    text(-1.2,0,['180',char(176)])
    % text(0,-1.05,['-270',char(176)])
end
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