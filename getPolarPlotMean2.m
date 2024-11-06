function getPolarPlotMean2(data,color)

if nargin<2
    color = [];
end
if isempty(color)
    color = parula(numel(data));
end

%create the mean direction plot
for i =1:numel(data)
    phaseValues = reshape(data{i},[numel(data{i}),1]);
    phaseValues(isnan(phaseValues)) = [];
    dataMean = circ_mean(phaseValues);
    quiver(0,0,cos(dataMean),sin(dataMean),'AutoScaleFactor',0.9,'LineWidth',3,'Color',color(i,:))
    hold on
end    
%
axis square
line([-1 1],[0 0],'Color',[0.811764705882353,0.811764705882353,0.811764705882353])
line([0 0],[-1 1],'Color',[0.811764705882353,0.811764705882353,0.811764705882353])
circle(0,0,1,[0,0,0]);
circle(0,0,0.5,[0.811764705882353,0.811764705882353,0.811764705882353]);
axis([-1,1,-1,1])
text(1.01,0,['0',char(176)])
text(-1.2,0,['180',char(176)])
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