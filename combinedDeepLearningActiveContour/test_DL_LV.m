% testing Deep Learning algorithms
clear all
clc
close all
addpath('functions');
% save('DL_ROI_params','stackedAEOptTheta','inputSize','hiddenSizeL1','hiddenSizeL2','outputSize','netconfig');
% save('DL_LV_params','stackedAEOptTheta','inputSize','hiddenSizeL1','hiddenSizeL2','outputSize','netconfig');
%%
% load images
load matFiles/validation_data; 

patient=15;
%min_sn=sum(slice_per_patient(1:patient-1))+1;
min_sn=1;
max_sn=sum(slice_per_patient(1:patient-1))+slice_per_patient(patient);
slice_num=min_sn:max_sn;

% get dicom info
dicom_path=['dcom/validation/ED',num2str(patient),'/images'];
para=get_dicominfo(dicom_path);


% show plots or not
disp_ena=0;

% load parameters of Deep Learning for ROI
Mroi=100;


% load Deep learning paramenters for LV segmentation
load DL_LV_params.mat;

for k=1:length(slice_num)
%for k=1:size(t_Iroi,3)
    
    C1=t_contours{k}; % read manual contours
    m_cnt=t_centers{k}; % read centers
    %subI=t_Iroi(:,:,k);
    display(['processing slice# ',num2str(slice_num(k))])
    subI=t_Iroi(:,:,slice_num(k));

    
    % run Deep Learning network to find the segmentation of LV
    t_yLV_h(:,:,k)=DLN(subI,stackedAEOptTheta,inputSize,hiddenSizeL1,hiddenSizeL2,outputSize,netconfig);
    
    % clean segmentations
    t_yLV_clean(:,:,k)=clean_segs(t_yLV_h(:,:,k));
    
    % dice metric
    [overlap, dice(k)]=DiceSimilarity2DImage(t_yLV_clean(:,:,k), t_yLV(:,:,k));
    
   
    % resize to original image size
    t_LVmask(:,:,k)=remap_mask(t_yLV_clean(:,:,k),t_centers{k},t_yROI(:,:,k));
    autoPoints=contourc(double(t_LVmask(:,:,k)), [0 0]);autoPoints=autoPoints(:,2:end)';

    % Perpendicular Distance
    manualPoints=t_contours{slice_num(k)};
    PD(k) = calc_dist(autoPoints,manualPoints,para);    

end

%% dice metric
avg_dice=mean(dice)
perp_dis=mean(PD)

%% display images and contours

if disp_ena==1
% get number of studies
num_studies=length(slice_per_patient);

for k=1:num_studies

    % get number of slices per patient
    spp=slice_per_patient(k);
  
    for k2=1:spp
        % get manual contour file name
        ind=k2+sum(slice_per_patient(1:k-1));
        figure(k)
        subplot(5,3,k2)
        imagesc(t_Iroi(:,:,ind));colormap(gray);
        hold on; contour(t_yLV_clean(:,:,ind),[0 0],'g','LineWidth',2); 
        contour(t_yLV(:,:,ind),[0 0],'r','LineWidth',2); 
        
        figure(num_studies+k)
        subplot(5,3,k2)
        imagesc(t_I(:,:,ind));colormap(gray);
        hold on; contour(t_LVmask(:,:,ind),[0 0],'r','LineWidth',2); 
        plot(t_contours{ind}(:,1),t_contours{ind}(:,2),'g','LineWidth',2); 
        
    end     
end
end
%% save contours as text files
% output=save_contours(t_LVmask,t_LV_cont_names,slice_per_patient);

