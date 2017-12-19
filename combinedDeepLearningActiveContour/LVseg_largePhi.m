% two-dimensional segmenter  combined with Deep Learning
clear all   
close all
clc
addpath('functions')
%% load the image 
disp('Load MRI images');
load matFiles/validation_data; 
%%
slice_num=23;
I=t_I(:,:,slice_num);
figure(1)
maxplots=4;
subplot(1,maxplots,1)
imagesc(I);
colormap(gray)
title(['slice number=',num2str(slice_num)]);

% maximum iteration
max_its = 50;
intEweight=.5;
DLweight=0.2;
Dynamic_Window=1;

% region of interest size
Mroi=100;

%--- run DL to find ROI: LV location detection
% load parameters of Deep Learning for ROI
%load DL_ROI_params.mat;
%yROI=DLN(I,stackedAEOptTheta,inputSize,hiddenSizeL1,hiddenSizeL2,outputSize,netconfig);
%subplot(1,maxplots,2)
%magesc(yROI);
%subplot(1,maxplots,2);hold on; plot(m_cnt(1),m_cnt(2),'r+')

% cut ROI from the image
%[subI,m_cnt]=mask2subImage(I,yROI,Mroi);
subI=t_Iroi(:,:,slice_num);
m_cnt=t_centers{slice_num};
subplot(1,maxplots,3)
imagesc(subI);colormap(gray);hold on; plot(50,50,'r+')


% run DL-LV segmentation to get phi_0
load DLconfigure/DL_LV_params.mat;
mask_LV=DLN(subI,stackedAEOptTheta,inputSize,hiddenSizeL1,hiddenSizeL2,outputSize,netconfig);
mask_LV=clean_segs(mask_LV);
init_mask = remap_mask(mask_LV,m_cnt,I);

% run segmentation
figure(2)
[LV_seg,phi,m_cnt] = region_segLargePhi(I,subI,m_cnt,init_mask,max_its,intEweight,DLweight,Dynamic_Window,1);
hold on
%imshow(subI,'initialmagnification',200,'displayrange',[0 255]); hold on;
contour(init_mask, [0 0], 'y','LineWidth',2);
Ct=t_contours{slice_num};
plot(Ct(:,1),Ct(:,2),'r','LineWidth',2)
legend('auto','auto','initial','mamual');




