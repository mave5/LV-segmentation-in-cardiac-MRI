% read dcom images and collect all of x y z of data
clc 
close all
clear all
addpath('functions');
%imset='validation';
%imset='training';
imset='online';




%%
for patient=1:15
    %read direcotory of images
    folder=['dcom/',imset,'/ED',num2str(patient),'/images'];
    [~, xthickness(patient), ythickness(patient), zthickness(patient)] = gatherImages(folder);
end

patient_nums=1:15;
A= [patient_nums',xthickness',ythickness',zthickness'];

filename=['dcom/',imset,'/xyz_info'];
xlswrite(filename,A);