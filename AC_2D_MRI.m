% LV segmentation from 2D cardiac MRI
% M. R. Avendi, 2014-2015

clear all   
close all
clc
addpath('functions')
%% load the image 
disp('Load MRI images');
load ('matFiles/images1.mat','Iroi','yLV')

max_its = 200; % maximum iterations
intEweight=0.5; % weight of length energy 
ShapeWeight=0.0;  % weight of shape engery


slice_num=1; % slice number
subI=Iroi(:,:,slice_num);
ground_truth=yLV(:,:,slice_num);

% show image
showCurveAndPhi(subI,ground_truth);
legend('ground truth')
disp('Please draw an initial contour...');

% initialization
h = imfreehand(gca,'closed',false); 
init_mask=createMask(h);

% run segmentation
disp('sementation in progress ...');
[auto_seg1,phi] = ac_seg(subI,init_mask,max_its,intEweight,ShapeWeight,1);

% clean segmentation, remove islands and small contours
auto_seg2=clean_segs(auto_seg1);

% show automatic and manual segmentation
showCurveAndPhi(subI,ground_truth,auto_seg2);
legend('ground truth','','automatic')
disp('sementation completed.');