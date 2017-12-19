% read mat files and store them as a 3D matrix for further processing
clc 
close all
clear all
%%
imset='training';
%imset='validation';
%imset='online';

folder=['matFiles/',imset,'/ED/'];
D=dir([folder,'\*.mat']);
numpatients = length(D(not([D.isdir])));

t_I=[];
t_yROI=[];
t_Iroi=[];
t_yLV=[];
t_yMC=[];
t_contoursLV=[];
t_contoursMC=[];
t_centers=[];
t_LV_cont_names=[];
t_MC_cont_names=[];
for k=1:numpatients
    imgname=['matFiles/',imset,'/ED/ED_P',num2str(k)];
%    yname=['matFiles/',imset,'/yROI_ED',num2str(k)];
    load (imgname);
%    load (yname);
    t_I=cat(3,t_I,I);
    t_yROI=cat(3,t_yROI,yROI);
    t_Iroi=cat(3,t_Iroi,Iroi);
    
    t_yLV=cat(3,t_yLV,yLV);
    t_yMC=cat(3,t_yMC,yMC);
    
    t_contoursLV=[t_contoursLV,contours_LV(1:size(I,3))];
    t_contoursMC=[t_contoursMC,contours_MC(1:size(I,3))];
    
    t_centers=[t_centers,contour_center(1:size(I,3))];
    slice_per_patient(k)=size(I,3);
    t_LV_cont_names=[t_LV_cont_names;LV_contours_names];
    t_MC_cont_names=[t_MC_cont_names;MC_contours_names];
end

filename=['matFiles/',imset,'_dataED'];

save (filename,'t_I','t_yROI','t_Iroi','t_yLV','t_yMC','t_contoursLV','t_contoursMC','t_centers','numpatients',...
'slice_per_patient','t_LV_cont_names','t_MC_cont_names');
    
display('file has been saved')    