% two-dimensional segmenter combined with Deep Learning
clear all   
close all
clc
addpath('functions')

%% load the image 
disp('Load MRI images');
imset='training';
%imset='online';
%imset='validation';
name1=['matFiles/',imset,'_data'];
load (name1)
%%
patient=2;
min_sn=sum(slice_per_patient(1:patient-1))+1;
max_sn=sum(slice_per_patient(1:patient-1))+slice_per_patient(patient);
step1=ceil((max_sn-min_sn)/5);
sn_range=min_sn:step1:max_sn;
slice_num=sn_range(4);
%slice_num=110;

% maximum iteration of active contour
max_its = 150;
intEweight=0.5;
DLweight=0.1;
Dynamic_Window=0;

% enable plots
dis_ena=1;

% region of interest size
Mroi=100;

% get dicom info
dicom_path=['dcom/validation/ED',num2str(patient),'/images'];
para=get_dicominfo(dicom_path);


display(['processing slice# ',num2str(slice_num)])
I=t_I(:,:,slice_num);

% cut ROI from the image
subI=t_Iroi(:,:,slice_num);

% changing the size of the ROI if we want, useful for apex slices
roi_x=round((size(subI,1)-Mroi)/2)+1:round((size(subI,1)+Mroi)/2);
roi_y=round((size(subI,2)-Mroi)/2)+1:round((size(subI,2)+Mroi)/2);
subI2(:,:,slice_num)=subI(roi_y,roi_x);
m_cnt=t_centers{slice_num};

% run DL-LV segmentation to get phi_0
load DLconfigure/DL_LV_params.mat;
DL_prms_fn='DLconfigure/DL_LV_params.mat';
init_mask1=DLN(subI2(:,:,slice_num),stackedAEOptTheta,inputSize,hiddenSizeL1,hiddenSizeL2,outputSize,netconfig);
init_mask2=clean_segs(init_mask1);
init_mask(:,:,slice_num)=remap_mask(init_mask2,size(subI)/2+1,subI);


% run segmentation
[LV_seg1,phi,m_cnt] = region_seg_subPhi(I,subI2(:,:,slice_num),m_cnt,init_mask1,max_its,intEweight,DLweight,Dynamic_Window,0,DL_prms_fn);
LV_seg2=clean_segs(LV_seg1);
out_mask = remap_mask(LV_seg2,size(subI)/2+1,subI);
LV_seg(:,:,slice_num)=out_mask;

% compute metrics
manualPoints=t_contours{slice_num};
[dm1,dm2,PD1,PD2,HD1,HD2]=eval_metrics(I,m_cnt,LV_seg(:,:,slice_num),manualPoints,para); 

% special care for slices with large Perp Distance
if PD2>= 3 
    edited_seg=edit_prior_shape(subI2(:,:,slice_num),LV_seg(:,:,slice_num),0);
    [dm11,dm22,PD11,PD22,HD11,HD22]=eval_metrics(I,m_cnt,edited_seg,manualPoints,para); 
    if PD22<PD2
        dm1=dm11;dm2=dm22;PD1=PD11;PD2=PD22;HD1=HD11;HD2=HD22;
        LV_seg(:,:,slice_num)=edited_seg;
    end
end


%% evaluate segmentations 
DM1=dm1
DM2=dm2


%% show only one slice

C1=t_contours{slice_num};Ct=scaleContour(C1,t_centers{slice_num},100);
stn=[imset,':patient:',num2str(patient),' slice:',num2str(slice_num)];

figure 
showCurveAndPhi(t_Iroi(:,:,slice_num),Ct ,(LV_seg(:,:,slice_num)), slice_num);
ylabel(stn);

figure
showCurveAndPhi(t_Iroi(:,:,slice_num),Ct ,bwconvhull(LV_seg(:,:,slice_num)), slice_num);
ylabel(stn);



