% dimension reduction using PCA
clc
clear all
close all
addpath('functions')
%%
load matFiles/training_data; 

% total images
I=t_Iroi;
[x,y,z]=size(I);
I1=reshape(I,x*y,z);

% apply PCA
N=size(I1,1);
C=I1*I1'/N;
[U,S,V]=svd(C);

% reduce dimension
n=90;
I2=I1'*U(:,1:n);

%save pca_results
%% recover
I3=I2*U(:,1:n)';
I4=I3';
I5=reshape(I4,x,y,z);
imagesc(I5(:,:,1));colormap(gray)

% total ROI
B=t_yROI;
[x,y,z]=size(B);
B1=reshape(B,x*y,z);

p2=pca(B1');
B2=B1'*p2;

% recover
B3=B2*p2';
B4=B3';
B5=reshape(B4,x,y,z);
figure
imagesc(B5(:,:,2));colormap(gray)
imagesc(B(:,:,2));colormap(gray)
Br=imresize(B(:,:,2),40/256);
Brr=imresize(Br,256/40);
Brr=round(Brr);
imagesc(Brr);colormap(gray)


