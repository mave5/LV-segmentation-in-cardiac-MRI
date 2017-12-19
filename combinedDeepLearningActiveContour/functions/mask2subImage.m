% this function takes a mask, an image I, and a size M, and extracts an M*M
% sub-image I_sub centered at cnt.

function [I_sub,cnt]=mask2subImage(I,mask,M)
% inputs
% I : original image
% mask : mask should be the same size of the original image
% M : size of ROI

% output
% I_sub : sub image centered at cnt
% cnt : ceneter of the sub-image respect to the original image

% find the center of the mask
temp=regionprops(mask, 'Centroid');
cnt1=temp.Centroid;
%cnt=floor(cnt1);
cnt=round(cnt1);


% coordinates of the top and bottom corners
cnt_x=cnt(1);
cnt_y=cnt(2);

% top left corner
x1=cnt_x-M/2;
y1=cnt_y-M/2;

% bottom right corner
x4=cnt_x+M/2-1;
y4=cnt_y+M/2-1;

% extract the sub-image
I_sub=I(y1:y4,x1:x4);


end