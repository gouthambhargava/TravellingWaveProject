function [allPoints] = getRawOverlap(olBounds,olDirs1,olDirs2,req)
if nargin <4
    req = 0;
end
allPoints1 = [];
allPoints2 = [];
for i = 1:numel(olBounds(1,:))
    bounds1 = olBounds{1,i};
    bounds2 = olBounds{2,i};
    allOverlaps1 = [];
    allOverlaps2 = [];
    for j = 1:size(bounds1,2)
        olPoints = intersect(bounds1(1,j):bounds1(2,j),bounds2(1,j):bounds2(2,j));
        dir1 = olDirs1(i,olPoints);
        dir2 = olDirs2(i,olPoints);
        if req==1
            allOverlaps1 = cat(2,allOverlaps1,circ_mean(dir1));
            allOverlaps2 = cat(2,allOverlaps2,circ_mean(dir2));
        elseif req==2
            allOverlaps1 = cat(2,allOverlaps1,unique(dir1));
            allOverlaps2 = cat(2,allOverlaps2,unique(dir2));
        else
            allOverlaps1 = cat(2,allOverlaps1,dir1);
            allOverlaps2 = cat(2,allOverlaps2,dir2);
        end
    end
allPoints1 = cat(2,allPoints1,allOverlaps1);
allPoints2 = cat(2,allPoints2,allOverlaps2);
end
if ~req==2
    allPoints = cat(1,allPoints1,allPoints2);
else 
    allPoints = {allPoints1,allPoints2};
end
end
        
    