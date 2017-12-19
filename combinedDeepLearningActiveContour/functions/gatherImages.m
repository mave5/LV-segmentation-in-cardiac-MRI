function [threedarray, xthickness, ythickness, zthickness] = gatherImages(folder)
%GATHERIMAGES looks through a folder, gets all DICOM files and assembles
%them into a viewable 3d format.
currentfolder=pwd;
d = sortDirectory(folder); %Sort in ascending order of instance number
topimage = dicomread(cell2mat(d(1,:)));
metadata = dicominfo(cell2mat(d(1,:)));

[group1 element1] = dicomlookup('PixelSpacing');
[group2 element2] = dicomlookup('SliceThickness');
resolution = metadata.(dicomlookup(group1, element1));
xthickness = resolution(1); ythickness = resolution(2);
zthickness = metadata.(dicomlookup(group2, element2));

threedarray = zeros(size(topimage,1),size(topimage,2),size(d,1));
threedarray(:,:,1) = topimage;

for k2 = 2:size(d,1)
   threedarray(:,:,k2) = dicomread(cell2mat(d(k2,:)));
end
cd (currentfolder);


function d = sortDirectory(folder)
%SORTDIRECTORY sorts based on instance number

cd(folder);
%d = ls('*.dcm');
dfile = dir('*.dcm'); 
m = size(dfile,1);

[group, element] = dicomlookup('InstanceNumber');
sdata(m) = struct('imagename','','instance',0);

for k1 = 1:m
    metadata = dicominfo(dfile(k1,:).name);
    position = metadata.(dicomlookup(group, element));
    sdata(k1) = struct('imagename',dfile(k1,:).name,'instance',position);
end

[unused order] = sort([sdata(:).instance],'ascend');
sorted = sdata(order).';

for k1 = 1:m
    d(k1,:) = {sorted(k1).imagename};  % we use cell in case of different filename length.
end





