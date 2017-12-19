% display gray images
function disImgs(I,num)
% I : input images as a 2D/3D matrix
% number of images to displayed

if nargin==1
    num=1;
end

[x y z]=size(I);

npx=ceil(num/5);
npy=max([ceil(num/5),5]);

figure
for k=1:num
    subplot(npx,npy,k);
    imagesc(I(:,:,k));
    colormap(gray);
end


end