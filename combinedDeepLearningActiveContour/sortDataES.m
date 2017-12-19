% read mat files and store them as a 3D matrix for further processing
clc 
close all
clear all
%%
%imset='training';
imset='validation';
%imset='online';

folder=['matFiles/',imset,'/ES/'];
D=dir([folder,'\*.mat']);
numpatients = length(D(not([D.isdir])));

t_I=[];
t_yROI=[];
t_Iroi=[];
t_yLV=[];
t_contours=[];
t_centers=[];
t_LV_cont_names=[];
t_MC_cont_names=[];
for k=1:numpatients
    imgname=['matFiles/',imset,'/ES/ES_P',num2str(k)];
%    yname=['matFiles/',imset,'/yROI_ED',num2str(k)];
    load (imgname);
%    load (yname);
    t_I=cat(3,t_I,I);
    t_yROI=cat(3,t_yROI,yROI);
    t_Iroi=cat(3,t_Iroi,Iroi);
    t_yLV=cat(3,t_yLV,yLV);
    t_contours=[t_contours,contours(1:size(I,3))];
    t_centers=[t_centers,contour_center(1:size(I,3))];
    slice_per_patient(k)=size(I,3);
    t_LV_cont_names=[t_LV_cont_names;LV_contours_names];
end

filename=['matFiles/',imset,'_dataES'];

save (filename,'t_I','t_yROI','t_Iroi','t_yLV','t_contours','t_centers','numpatients',...
'slice_per_patient','t_LV_cont_names');
    
display('file has been saved')    