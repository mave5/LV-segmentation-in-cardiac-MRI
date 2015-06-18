
% show image, contours
function showCurveAndPhi(varargin)

I=varargin{1};
max_range=min(255,max(I(:)));
imshow(I,'initialmagnification',200,'displayrange',[0 max_range]); hold on;
colsty1=['g';'r';'b';'k';'y'];
for k=1:nargin-1
  
    B1=varargin{k+1};
    if size(B1,1)==2 || size(B1,2)==2
       plot(B1(:,1),B1(:,2),colsty1(k),'LineWidth',3);
    else
       contour(B1, [0 0], colsty1(k),'LineWidth',3);
       contour(B1, [0 0], 'k','LineWidth',1);
    end
end

end
