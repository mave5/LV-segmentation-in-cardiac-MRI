% testing Deep Learning ROI
clear all
clc
close all
addpath('functions');
%save('DL_ROI_params','stackedAEOptTheta','inputSize','hiddenSizeL1','hiddenSizeL2','outputSize','netconfig');
%%

% load parameters of Deep Learning for ROI
Mroi=100;

% load images
load matFiles/validation_data; 

% load DL ROI
load DLconfigure/DL_ROI_params.mat;

% choose a slice
nsl=28;
for slice_num=1:nsl

I=t_I(:,:,slice_num);
%tt1=randi(10)*0;
%I = imtranslate(I, [tt1 tt1]);
C1=t_contours{slice_num};

% show image
%subplot(6,5,slice_num)
%imagesc(I); colormap(gray);
%hold on
%plot(C1(:,1),C1(:,2),'r');

% run DL to find ROI: LV location detection
t_yROI_h(:,:,slice_num)=DLN(I,stackedAEOptTheta,inputSize,hiddenSizeL1,hiddenSizeL2,outputSize,netconfig);
%subplot(5,6,slice_num)
%imagesc(t_yROI_h(:,:,slice_num));colormap(gray);

% cut ROI from the image
[t_Iroi_h(:,:,slice_num),m_cnt(slice_num,:)]=mask2subImage(I,t_yROI_h(:,:,slice_num),Mroi);
subplot(5,6,slice_num)
imagesc(t_Iroi_h(:,:,slice_num));colormap(gray);

% find center from manual contours
Cx_LV=C1(:,1);Cy_LV=C1(:,2);
[junk,cnt_man(slice_num,1),cnt_man(slice_num,2)]=polycenter(Cx_LV,Cy_LV);

end

save ('test_ROI.mat','t_Iroi_h','m_cnt');

% compute dice metric
t_dice=0;
for k=1:size(t_yROI_h,3)
    img1=t_yROI_h(:,:,k);
    img2=t_yROI(:,:,k);
    [overlapI,dice(k)]=DiceSimilarity2DImage(img1, img2);
end
average_dice=mean(dice)

%%
% slice_num=1;
% subI=t_Iroi_h(:,:,slice_num);
% otsuMask=otsu(subI);
% imagesc(otsuMask);colormap(gray)
%otsuMask_clean=clean_segs(otsuMask);
%imagesc(otsuMask_clean);colormap(gray)

% otsuMask=otsuMask>1;
% subplot(1,2,1)
% imagesc(otsuMask);colormap(gray)
% CC = bwconncomp(otsuMask);
% L = labelmatrix(CC);
% L1=L==5;
% subplot(1,2,2)
% imagesc(L1);colormap(gray);
% 
% % find the center of the mask
% temp=regionprops(CC, 'Centroid');
% cnt1=temp.Centroid;
% cnt=round(cnt1);
