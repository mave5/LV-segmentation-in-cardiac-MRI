% this function takes the image, and optimized parameters of the deep
% learning network and outputs a mask

function y=DLN(I,parameters,inputSize,hiddenSizeL1,hiddenSizeL2,outputSize,netconfig)
% I : image to be masked
% parameters: the learned/optimized parameters 
% inputSize : the visible size
% hiddenSizeL1 and L2: the size for layer 1 and 2
% outputSize
% netconfig

patchsize=sqrt(inputSize);
test_input=sampleIMAGES(I,patchsize);

stackedAEOptTheta=parameters;

[pred_y] = stackedAEPredict(stackedAEOptTheta, inputSize, hiddenSizeL1, ...
                          outputSize, netconfig, test_input);

% mask output
y=reshape(pred_y,patchsize,patchsize,[]);

scale=size(I,1)/patchsize;
y=imresize(y,scale);

end