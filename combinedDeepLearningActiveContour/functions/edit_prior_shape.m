
% prior shape might be inaccurate for bottom slices
% to get better prior shape, we use intersection between otsu and prior
% shape

function CH=edit_prior_shape(I,prior,disp_ena)

prior=logical(bwconvhull(prior));

I1=I.*prior;

% R1=40:60;
% I2=I1(R1,R1);
% imshow(I2,[0 255])
% ot1=otsu(I2)>1;
% imshow(ot1)
% otsu thresholding
ot1=otsu(I1)>1;

% intersection of otsu and prior shape
com1=logical(ot1).*prior;

% find the region at center
cnt=round(size(I)/2);
com2=regiongrowing(com1,cnt(1),cnt(2));

if sum(com2(:))<5
    com2=prior;
end
%BW=clean_segs(com1);
BW=com2;

% covenx hull
CH = bwconvhull(BW);

if disp_ena==1
    imshow(I,'initialmagnification',200,'displayrange',[0 255]); hold on;
    contour(prior, [0 0], 'g','LineWidth',2);
    contour(BW, [0 0], 'r','LineWidth',2);
    contour(CH, [0 0], 'b','LineWidth',2);
    legend('original','edited','Convex Hull');
end

end

