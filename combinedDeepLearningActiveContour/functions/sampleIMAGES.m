function patches = sampleIMAGES(IMAGES,patchsize,norm_ena)
% sampleIMAGES
% inputs
% numpatches: for example 1E4
% patchsize% for example 8*8
% IMAGES: images in a 3D matrix
% output
% patches % a vector of randomly chosen patches
if nargin==2
    norm_ena=1;
end
visibleSize = patchsize*patchsize;   % number of input units 

% get size and number of images
[xn, yn, zn]=size(IMAGES);

scale=patchsize/xn;
imgs = imresize(IMAGES, scale);
patches=(reshape(imgs,visibleSize,zn));

%% ---------------------------------------------------------------
% For the autoencoder to work well we need to normalize the data
% Specifically, since the output of the network is bounded between [0,1]
% (due to the sigmoid activation function), we have to make sure 
% the range of pixel values is also bounded between [0,1]
if norm_ena==1
patches = normalizeData(patches);
end
end

%% ---------------------------------------------------------------
function patches = normalizeData(patches)

% Squash data to [0.1, 0.9] since we use sigmoid as the activation
% function in the output layer

% Remove DC (mean of images). 
patches = bsxfun(@minus, patches, mean(patches));

% Truncate to +/-3 standard deviations and scale to -1 to 1
pstd = 3 * std(patches(:));
patches = max(min(patches, pstd), -pstd) / pstd;

% Rescale from [-1,1] to [0.1,0.9]
patches = (patches + 1) * 0.4 + 0.1;

end
