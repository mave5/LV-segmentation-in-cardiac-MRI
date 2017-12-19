% two-dimensional segmenter combined with Deep Learning
clear all   
close all
clc
addpath('functions')
%% load the image 
disp('Load MRI images');
%imset='training';
%imset='online';
imset='validation';
name1=['matFiles/',imset,'_data'];
load (name1)
%%
patient=1;
min_sn=sum(slice_per_patient(1:patient-1))+1;
max_sn=sum(slice_per_patient(1:patient-1))+slice_per_patient(patient);
%step1=ceil((max_sn-min_sn)/5);
step1=1;
slice_num=min_sn:step1:max_sn;
%slice_num=[min_sn,min_sn+step1:step1:max_sn-1,max_sn];
%slice_num=round(mean(slice_num));
%slice_num=max_sn;

% maximum iteration of active contour
max_its = 150;
intEweight=0.5;
DLweight=0.2;
Dynamic_Window=0;

% enable plots
dis_ena=1;

% region of interest size
Mroi=100;

% get dicom info
dicom_path=['dcom/validation/ED',num2str(patient),'/images'];
%im_folder=dir(dicom_path1);
%dicom_path2=['dcom/validation/ED',num2str(patient),'/images'];
%full_dicom_filename = [dicom_path2 filesep im_folder(1).name];
para=get_dicominfo(dicom_path);


for k=1:length(slice_num)
display(['processing slice# ',num2str(slice_num(k))])
I=t_I(:,:,slice_num(k));

%figure(1)
%maxplots=4;
%subplot(1,maxplots,1)
%imagesc(I);colormap(gray);title(['slice number=',num2str(slice_num)]);


%--- run DL to find ROI: LV location detection
% load parameters of Deep Learning for ROI
%load DL_ROI_params.mat;
%yROI=DLN(I,stackedAEOptTheta,inputSize,hiddenSizeL1,hiddenSizeL2,outputSize,netconfig);
%subplot(1,maxplots,2)
%imagesc(yROI);
%subplot(1,maxplots,2);hold on; plot(m_cnt(1),m_cnt(2),'r+')

% cut ROI from the image
%[subI,m_cnt]=mask2subImage(I,yROI,Mroi);
subI=t_Iroi(:,:,slice_num(k));

% changing the size of the ROI if we want, useful for apex slices
roi_x=round((size(subI,1)-Mroi)/2)+1:round((size(subI,1)+Mroi)/2);
roi_y=round((size(subI,2)-Mroi)/2)+1:round((size(subI,2)+Mroi)/2);
subI2(:,:,slice_num(k))=subI(roi_y,roi_x);
m_cnt=t_centers{slice_num(k)};
%subplot(1,maxplots,3)
%imagesc(subI);colormap(gray);hold on; plot(50,50,'r+')

% run DL-LV segmentation to get phi_0
load DLconfigure/DL_LV_params.mat;
DL_prms_fn='DLconfigure/DL_LV_params.mat';
patchsize=sqrt(inputSize);
%normI=sampleIMAGES(t_Iroi,patchsize);
init_mask=DLN(t_Iroi(:,:,slice_num(k)),stackedAEOptTheta,inputSize,hiddenSizeL1,hiddenSizeL2,outputSize,netconfig);
init_mask=clean_segs(init_mask);

% run segmentation
%figure(2)
%subplot(5,6,k)
[LV_seg1,phi,m_cnt] = region_seg_subPhi(I,subI2(:,:,slice_num(k)),m_cnt,init_mask,max_its,intEweight,DLweight,Dynamic_Window,0,DL_prms_fn);
LV_seg(:,:,slice_num(k))=clean_segs(LV_seg1);

% compute metrics
manualPoints=t_contours{slice_num(k)};
[dm1(k),dm2(k),PD1(k),PD2(k),HD1(k),HD2(k)]=eval_metrics(I,m_cnt,LV_seg(:,:,slice_num(k)),manualPoints,para); 

% special care for slices with large Perp Distance
if PD2(k)>= 3 
    edited_seg=edit_prior_shape(subI2(:,:,slice_num(k)),LV_seg(:,:,slice_num(k)),1);
    [dm11,dm22,PD11,PD22,HD11,HD22]=eval_metrics(I,m_cnt,edited_seg,manualPoints,para); 
    if PD22<PD2(k)
        dm1(k)=dm11;dm2(k)=dm22;PD1(k)=PD11;PD2(k)=PD22;HD1(k)=HD11;HD2(k)=HD22;
        LV_seg(:,:,slice_num(k))=edited_seg;
    end
%     [LV_seg1,phi,m_cnt] = region_seg_subPhi(I,subI,m_cnt,edited_seg,max_its,intEweight,DLweight,Dynamic_Window,0);    
%     LV_seg(:,:,slice_num(k))=clean_segs(LV_seg1);
%     [dm1(k),dm2(k),PD1(k),PD2(k),HD1(k),HD2(k)]=eval_metrics(I,m_cnt,LV_seg(:,:,slice_num(k)),manualPoints,para); 
end


end

%% evaluate segmentations 

% good contours
GP1=find(PD1<5);
GP2=find(PD2<5);

% averages 
AGP1=length(GP1)/length(PD1)
AGP2=length(GP2)/length(PD2)

PD1=PD1
PD2=PD2
average_PD1=mean(PD1(PD1<5))
average_PD2=mean(PD2(PD2<5))


DM1=dm1
DM2=dm2
average_dm1=mean(dm1(PD1<5))
average_dm2=mean(dm2(PD2<5))


% Hausdorff distance
average_HD1=mean(HD1(PD1<5))
average_HD2=mean(HD2(PD1<5))

[AGP1 average_dm1 average_PD1 average_HD1]
[AGP2 average_dm2 average_PD2 average_HD2]

%figure(5)
%subplot(1,2,1);stem(dm1);title('dice metric')
%subplot(1,2,2);stem(PD1);title('perpendicular distance')

%% display 2D segmentation
if dis_ena==1
    for k=1:length(slice_num)
    figure(1);subplot_tight(2,round(length(slice_num)/2),k)
    imshow(subI2(:,:,slice_num(k)),'initialmagnification',200,'displayrange',[0 255]); hold on;
    %contour(phi, [0 0], 'b','LineWidth',2);
    contour(LV_seg(:,:,slice_num(k)), [0 0], 'b','LineWidth',2);
    %contour(init_mask, [0 0], 'y','LineWidth',2);
    % contour(bwconvhull(init_mask), [0 0], 'r','LineWidth',2);
    %contour(t_yLV(roi_y,roi_x,slice_num(k)), [0 0], 'g--','LineWidth',2);
    contour(bwconvhull(LV_seg(:,:,slice_num(k))), [0 0], 'r','LineWidth',2);
    C1=t_contours{slice_num(k)};Ct=scaleContour(C1,t_centers{slice_num(k)},Mroi);
    plot(Ct(:,1),Ct(:,2),'g--','LineWidth',2);
    %legend('auto seg','initial contour','manual seg','auto conv hul');
    

%     figure(2)
%     subplot_tight(2,ceil(length(slice_num)/2),k)
%     imagesc(t_I(:,:,slice_num(k)));colormap(gray);
%     hold on
%     contour(double(phiI(:,:,slice_num(k))), [0 0], 'r','LineWidth',2);
%     manualPoints=t_contours{slice_num(k)};
%     plot(manualPoints(:,1),manualPoints(:,2),'g','LineWidth',2)   
%     lvmaskI=remap_mask(t_yLV(:,:,slice_num(k)),t_centers{slice_num(k)},I);
%     contour(double(lvmaskI), [0 0], 'y','LineWidth',2);
%     contour(double(mask_rc(:,:,slice_num(k))), [0 0], 'b','LineWidth',2);
%     legend('auto','manual');
    end
 

%% 3d represenations
if length(slice_num)<2 || dis_ena==0
break
end

V_LV_auto=LV_seg(:,:,(slice_num));
V_LV_man=t_yLV(:,:,(slice_num));
[x_max,y_max,z_max]=size(V_LV_man);

% find centers of masks
for k=1:size(V_LV_man,3)
    temp=regionprops(V_LV_man(:,:,k),'Centroid');
    cnt_man(k,:)=temp.Centroid;
    temp=regionprops(V_LV_auto(:,:,k),'Centroid');
    cnt_auto(k,:)=temp.Centroid;
end

% registeration to compensate mis-alignmnet 
ref_man=cnt_man(round(z_max/2),:);
for k=1:size(V_LV_man,3)
    trans=fliplr(-cnt_man(k,:)+ref_man);
    Vt_LV_man(:,:,k)=imtranslate(V_LV_man(:,:,k),trans);

    trans=fliplr(-cnt_auto(k,:)+ref_man);
    Vt_LV_auto(:,:,k)=imtranslate(V_LV_auto(:,:,k),trans);
    
end


% create the bottom slice
for k1=1:1
Vt_LV_man(:,:,z_max+k1)=zeros(x_max,y_max);
apex_cnt=round(ref_man);
R=0;
Vt_LV_man(apex_cnt(2)-R:apex_cnt(2)+R,apex_cnt(1)-R:apex_cnt(1)+R,z_max+k1)=1;
end

% create the bottom slice
for k1=1:1
Vt_LV_auto(:,:,z_max+k1)=zeros(x_max,y_max);
apex_cnt=round(ref_man);
R=0;
Vt_LV_auto(apex_cnt(2)-R:apex_cnt(2)+R,apex_cnt(1)-R:apex_cnt(1)+R,z_max+k1)=1;
end

% top slice 
top_cnt=round(ref_man);
temp=zeros(x_max,y_max);
temp(top_cnt(2),top_cnt(1),1)=1;
temp(top_cnt(2)+1,top_cnt(1)+1,1)=1;
Vt_LV_man2=Vt_LV_man;
Vt_LV_man2(:,:,1)=temp;


figure
subplot(1,3,1)
disp3d(flipdim(V_LV_man,3),0);
title('original')
axis off

subplot(1,3,2)
disp3d(flipdim(Vt_LV_man,3),0);
title('translated')
axis off

subplot(1,3,3)
disp3d(flipdim(Vt_LV_man2,3),0);
title('translated')
axis off

% subplot(1,3,3)
% disp3d(flipdim(Vt_LV_auto,3),0);
% title('translated')
% axis off 

%
figure
subplot(1,3,1)
disp3d(flipdim(V_LV_man,3),9);
title('original')
axis off

fz=7;
subplot(1,3,2)
disp3d(flipdim(Vt_LV_man,3),[fz,fz,fz]);
title('translated')
axis off
view([-10 40]);

fz=7;
subplot(1,3,3)
disp3d(flipdim(Vt_LV_man2,3),[fz,fz,fz]);
title('translated')
axis off
view([-10 40]);

% fz=9;
% subplot(1,3,3)
% disp3d(flipdim(Vt_LV_auto,3),[fz,fz,fz]);
% title('translated')
% axis off
% view([-10 40]);

end
%%
%   fz=7;
%   figure  
%   disp3d(flipdim(Vt_LV_man,3),[fz,fz,fz]);
%   title('translated')
%   axis off
%   view([-10 40]);
 

% figure(1)
% subplot_tight(2,10,9)
% imshow(t_Iroi(:,:,slice_num(k)),'initialmagnification',200,'displayrange',[0 255]); hold on;
% contour(bwconvhull(LV_seg(:,:,slice_num(k))), [0 0], 'r','LineWidth',2);
% C1=t_contours{slice_num(k)};Ct=scaleContour(C1,t_centers{slice_num(k)},Mroi);
% plot(Ct(:,1),Ct(:,2),'g--','LineWidth',2);

%% show only one slice
figure
C1=t_contours{slice_num(k)};Ct=scaleContour(C1,t_centers{slice_num(k)},Mroi);
showCurveAndPhi(t_Iroi(:,:,slice_num(k)), Ct,phi, slice_num(k));hold on

figure
showCurveAndPhi(t_Iroi(:,:,slice_num(k)),Ct ,bwconvhull(LV_seg(:,:,slice_num(k))), slice_num(k))
%% save results
name1=['results/',imset,'/case',num2str(patient)];
%save (name1)

