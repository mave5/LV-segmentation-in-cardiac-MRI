% read manual segmentation and generate a 3D plot
clc 
close all
clear all
addpath('functions');
%%
imset='training';
patient=1;
%% read direcotory of images
folder=['dcom/',imset,'/ED',num2str(patient),'/images'];
[IMAGES, xthickness, ythickness, zthickness] = gatherImages(folder);

 [x_max,y_max,z_max]=size(IMAGES);
% for k=1:z_max
%     subplot(3,5,k)
%     imagesc(IMAGES(:,:,k));
%     colormap(gray);
% end

% adjust the resolution so that it falls in range [0-max_res]
max_res=3;
norm_res = floor(max(IMAGES(:))/(max_res)) ; 
Im = (IMAGES./norm_res);

figure
for k1=1:size(Im,3)
I1=Im(:,:,k1);   
subplot(3,5,k1);
imshow(I1);
end  
%% read directory of contours
current_dir=pwd;
contour_dir=['dcom/',imset,'/ED',num2str(patient),'/contours/LV'];
cd(contour_dir);
dfile = dir('*.txt'); 
numFiles = size(dfile,1);

% load contours
for k=1:numFiles
    contours{k}=load (dfile(k,:).name);
end
cd(current_dir);

% display images and contours
for k=1:z_max
    figure (1)
    subplot(3,5,k)
    % show image
    I1=Im(:,:,k);
    imagesc(I1);
    colormap(gray);
    hold on

    % read contour
    C1=contours{k};   
    x1=floor(C1(:,1));
    y1=floor(C1(:,2));
    plot(x1,y1,'r');
    
    % creat segmentation mask
    figure (2)
    subplot(3,5,k)
    % note that row==y and column==x
    LV_seg(:,:,k)=accumarray([y1 x1],1,[y_max x_max]);    
    imagesc(LV_seg(:,:,k));
    colormap(gray)
    title('segmnetation mask')
end

%% find edges of segmentation contours
for k=1:size(Im,3)
    figure(3);
    subplot(3,5,k);
    %Find edges in intensity image
    edge_seg = edge(LV_seg(:,:,k)); 
    imshow(Im(:,:,k) -2*edge_seg); 
    title('edges')
end
%% show three dimensional patch
    figure
    phi = LV_seg(:,:,end:-1:1);
    p3 = patch(isosurface(phi, 0),  'FaceColor','red', 'EdgeColor','none'); 
    isonormals(phi,p3)
    zlabel('z','FontSize',10,'FontWeight','bold'); 
    box on ;  alpha(.6)
    xrange = []; 
    yrange = [];

    axis([1 size(phi,2) 1 size(phi,1) 1 size(phi,3)]) 
    set(gca,'YTick',yrange,'fontsize',8)
    set(gca,'XTick',xrange,'fontsize',8)
    set(gca,'ZTick',round(10:size(phi,3)/5:size(phi,3)), 'fontsize',8)
    
    hold on
    colormap(gray); 
    dimz=size(Im,3);
    for k1=1:3:dimz
    
    Img2D_CrossSection =  Im(:,:,k1); 
    xImage = [1 size(phi,2); 1 size(phi,2)]; %# The x data for the image corners
    yImage = [1 1 ; size(phi,1), size(phi,1)]; %# The y data for the image corners
    zImage = k1 * ones(2,2);   %# The z data for the image corners
    surf(xImage,yImage,zImage,...    %# Plot the surface
     'CData',Img2D_CrossSection ,...
     'FaceColor','texturemap'); 
    end    
    % viewpoint spec: view(az,el)
    view(-30,13)

axis tight
camlight
camlight(-80,-10)
lighting gouraud
