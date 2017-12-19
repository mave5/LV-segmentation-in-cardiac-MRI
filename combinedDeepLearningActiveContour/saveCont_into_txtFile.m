% this program aims to save automatic contours into text files
% to be later used for evaluation and statistical analysis
clear all   
close all
clc
addpath('functions')
%imset='training';
%imset='online';
imset='validation';

%% ED contours
patient=15;

name1=['results/',imset,'/ED/case',num2str(patient),'.mat'];
load (name1,'LV_seg','Mroi','slice_num','slice_per_patient','t_I',...
                    't_centers','t_contours','t_LV_cont_names');

% convert LV mask to original image size
for k=1:length(slice_num)
display(['processing slice# ',num2str(slice_num(k))])
I=t_I(:,:,slice_num(k));

% center of sub image
cnt_xy=t_centers{slice_num(k)};

% resize to original size
LV_seg_Isize(:,:,slice_num(k)) = remap_mask(LV_seg(:,:,slice_num(k)),cnt_xy,I);

% save into files
fname1=['results/',imset,'/contours-auto/case',num2str(patient),'/'];
if exist(fname1,'dir')==0
    mkdir(fname1);
end
endo_cn=t_LV_cont_names(slice_num(k)).name;
endo_cnm = strrep(endo_cn, 'manual', 'auto');

name2=[fname1,endo_cnm];
save_contours2(LV_seg_Isize(:,:,slice_num(k)),name2);

% display automatic and manual contorus
subplot(3,ceil(length(slice_num)/3),k)
Cm=t_contours{slice_num(k)};
showCurveAndPhi(I,Cm,LV_seg_Isize(:,:,slice_num(k)));

end
%%
% ES contours
name1=['results/',imset,'/ES/case',num2str(patient),'.mat'];
load (name1,'LV_seg','Mroi','slice_num','slice_per_patient','t_I',...
                    't_centers','t_contours','t_LV_cont_names');

figure
% convert LV mask to original image size
for k=1:length(slice_num)
display(['processing slice# ',num2str(slice_num(k))])
I=t_I(:,:,slice_num(k));

% center of sub image
cnt_xy=t_centers{slice_num(k)};

% resize to original size
LV_seg_Isize(:,:,slice_num(k)) = remap_mask(LV_seg(:,:,slice_num(k)),cnt_xy,I);

% save into files
fname1=['results/',imset,'/contours-auto/case',num2str(patient),'/'];
endo_cn=t_LV_cont_names(slice_num(k)).name
endo_cnm = strrep(endo_cn, 'manual', 'auto');

name2=[fname1,endo_cnm];
save_contours2(LV_seg_Isize(:,:,slice_num(k)),name2);

% display automatic and manual contorus
subplot(3,ceil(length(slice_num)/3),k)
Cm=t_contours{slice_num(k)};
showCurveAndPhi(I,Cm,LV_seg_Isize(:,:,slice_num(k)));

end
