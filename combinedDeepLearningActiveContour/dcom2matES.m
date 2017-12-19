% read dcom images and contours and save them into mat files
clc 
close all
clear all
addpath('functions');
imset='validation';
%imset='training';
%imset='online';
%%
patient=4;

% region of interest
Mroi=100;
%% read direcotory of images
folder=['dcom/',imset,'/ES',num2str(patient),'/images'];
%folder=uigetdir;
[IMAGES, xthickness, ythickness, zthickness] = gatherImages(folder);

[x_max,y_max,z_max]=size(IMAGES);
scale=Mroi/x_max;
 for k=1:z_max
 subplot(3,5,k)
 imagesc(IMAGES(:,:,k));colormap(gray);
 end

%% read directory of contours
current_dir=pwd;
%contour_dir_MC=['dcom/',imset,'/ED',num2str(patient),'/contours/MC/'];
contour_dir_LV=['dcom/',imset,'/ES',num2str(patient),'/contours/LV/'];

cd(contour_dir_LV);
LV_contours_names = dir('*.txt'); 
numFiles = size(LV_contours_names,1);

% load contours of LV
for k=1:numFiles
    contours{k}=load (LV_contours_names(k,:).name);
end
cd(current_dir);


nLarge=z_max;
yROI=zeros(x_max,y_max,nLarge);
Iroi=zeros(Mroi,Mroi,nLarge);
for k=1:z_max

    % show images
    figure (1)
    subplot(3,5,k)
    I1=IMAGES(:,:,k);
    imagesc(I1);
    colormap(gray);
    title(['image',num2str(k)])
    hold on
   
    % read contour and plot it
    C_LV=contours{k};
    Cx_LV=C_LV(:,1);Cy_LV=C_LV(:,2);
    
    % find the contour center and plot it
    [junk,xcnt,ycnt]=polycenter(Cx_LV,Cy_LV);
    plot(Cx_LV,Cy_LV,'r','LineWidth',2);   
    plot(xcnt,ycnt,'r*','markerSize',12);

    % creat segmentation mask and display it
    figure(4)
    subplot(3,5,k)
    LV_seg(:,:,k)= roipoly(I1,Cx_LV, Cy_LV);    
    imagesc(LV_seg(:,:,k));
    colormap(gray)
    title(['mask',num2str(k)])
    hold on
    plot(Cx_LV,Cy_LV,'b','LineWidth',2)
    plot(xcnt,ycnt,'b*','markerSize',12); 
   
    x_cnt=round(xcnt);
    y_cnt=round(ycnt);
    contour_center{k}=[x_cnt,y_cnt];
    
    % define a rectangle centered at contour
    x_roi=x_cnt-Mroi/2:x_cnt+Mroi/2-1;
    y_roi=y_cnt-Mroi/2:y_cnt+Mroi/2-1;
    xroi=[x_cnt-Mroi/2,x_cnt+Mroi/2,x_cnt+Mroi/2,x_cnt-Mroi/2,x_cnt-Mroi/2];
    yroi=[y_cnt-Mroi/2,y_cnt-Mroi/2,y_cnt+Mroi/2,y_cnt+Mroi/2,y_cnt-Mroi/2];

    % create ROI mask: this will be used for training of DL-ROI
    figure(1)
    subplot(3,5,k)
    yROI(:,:,k)=poly2mask(xroi,yroi,x_max,y_max);
    contour(yROI(:,:,k),[0 0],'r')
    plot(xcnt,ycnt,'b*','markerSize',12);
    
    % find ROI in the mask: this will be used as training data for DL-LV
    figure(5)
    subplot(3,5,k)
    yLV(:,:,k)=LV_seg(y_roi,x_roi,k);
    imshow(yLV(:,:,k));
    title('yLV')
    hold on
    plot(50,50,'r*','markerSize',12);


    % ROI in the image: this will be used as input training data for DL-LV 
    figure(3)
    subplot(3,5,k)
    Iroi(:,:,k)=I1(y_roi,x_roi);
    imagesc(Iroi(:,:,k));
    colormap(gray);
    hold on
    plot(10,50,'r*','markerSize',12);
    contour(yLV(:,:,k),[0 0],'r','LineWidth',2)
    %plot(Cx_LV-x_cnt+Mroi/2+1,Cy_LV-y_cnt+Mroi/2+1,'r')
    
    % make sure that mask is ok
    %figure(6)
    %subplot(3,4,k)
    %I_LV(:,:,k)=LV_mask(:,:,k).*Iroi(:,:,k);
    %imagesc(I_LV(:,:,k));
    %title('LV')
    %colormap(gray);
    %hold on
    %plot(50,50,'r*','markerSize',12);

    % convert mask to poly for ROI
    %figure(6)
    %contour(LV_mask(:,:,k));
    %plot(xylv(:,1),xylv(:,2),'b')
    
end
%% store images and region of interest on disk
% P# stands for patient numebr
% ED stands for End Diastole
% yROI : this is a mask of ROI which is used as the output for traning the 
% DL ROI
% Iroi : this is part of the image that interests us and it will be used as
% the input of DL-LV segmentation
% yLV : this is a mask which is used as the output for training the DL-LV
% segmentation
I=IMAGES(:,:,1:nLarge);
filename=['matFiles/',imset,'/ES/ES_P',num2str(patient)];
save (filename, 'I', 'yROI','Iroi','yLV','contours','patient','xthickness', 'ythickness' ,'zthickness','contour_center','LV_contours_names');
disp('results saved')

% make sure sizes are matched
sizes=[size(I,3),size(Iroi,3),size(yLV,3),size(yROI,3)]
%%
%figure
%h1=disp3d(flipdim(LV_seg,3),'red',7);
%hold on
%h2=disp3d(flipdim(MC_seg,3),'g',7);
%alpha(.6)

%colormap(gray); 
%z1 = ceil(1);
%z1=z_max;
%I1 =  IMAGES(:,:,z1); 
%xImage = [1 y_max; 1 y_max]; %  The x data for the image corners
%yImage = [1 1 ; x_max, x_max]; %  The y data for the image corners
%zImage = (z_max-z1+1) * ones(2,2);   % The z data for the image corners

%surf(xImage,yImage,zImage,...    %  Plot the surface
%     'CData',I1 ,...
%     'FaceColor','texturemap'); 


