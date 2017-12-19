
% scale a given contour based on a center and widow size
function   Ct=scaleContour(C,cnt,M)
% C : input contour
% cnt : center
% M   : window size


Ctx=C(:,1)-cnt(1)+M/2+.5;
Cty=C(:,2)-cnt(2)+M/2+.5;

Ct=[Ctx,Cty];


end
