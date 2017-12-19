%% automatic detection of heart
clc;
clear all;
close all;
addpath('functions');
%% STEP 0: 
% parameters
patchsize = 32;
visibleSize = patchsize*patchsize;   % number of input units 
hiddenSizeL1 = 100;     % number of hidden units 
hiddenSizeL2=100;
sparsityParam1 = 0.01;   % desired average activation of the hidden units.
lambda = 3e-3;       % weight decay parameter       
beta = 3;            % weight of sparsity penalty term       
outputSize=visibleSize; % number of output units
%%======================================================================
%% STEP 1: laod training inputs and labels from mat files
load matFiles/training_data; 
train_input=sampleIMAGES(t_I,patchsize);
train_labels=sampleIMAGES(t_yROI,patchsize);
%% train sparse Auto Encoder 1
%  Randomly initialize the parameters
saeTheta1 = initializeParameters(hiddenSizeL1, visibleSize);

%  Use minFunc to minimize the function
addpath minFunc/
options.Method = 'lbfgs'; % Here, we use L-BFGS to optimize our cost
                          % function. Generally, for minFunc to work, you
                          % need a function pointer with two outputs: the
                          % function value and the gradient. In our problem,
                          % sparseAutoencoderCost.m satisfies this.
options.maxIter = 400;	  % Maximum number of iterations of L-BFGS to run 
options.display = 'on';

[sae1OptTheta, cost] = minFunc( @(p) sparseAutoencoderCost(p, ...
                                   visibleSize, hiddenSizeL1, ...
                                   lambda, sparsityParam1, ...
                                   beta, train_input), ...
                                   saeTheta1, options);
%% STEP 5: Visualization of AE1
W1 = reshape(sae1OptTheta(1:hiddenSizeL1*visibleSize), hiddenSizeL1, visibleSize);
display_network(W1', 12); 
%% compute activations from layer 1
[sae1Features] = feedForwardAutoencoder(sae1OptTheta, hiddenSizeL1, ...
                                        visibleSize, train_input);

%% train sparse Auto Encoder 2                                   
%  Randomly initialize the parameters
sae2Theta = initializeParameters(hiddenSizeL2, hiddenSizeL1);
sparsityParam2=1e-1;
lambda2 = 3e-3;       % weight decay parameter       
beta2 = 3;            % weight of sparsity penalty term       
[sae2OptTheta, costL2] = minFunc( @(p) sparseAutoencoderCost(p, ...
                                  hiddenSizeL1, hiddenSizeL2, ...
                                  lambda2, sparsityParam2, ...
                                  beta2, sae1Features), ...
                                  sae2Theta, options);
                              
W2 = reshape(sae2OptTheta(1:hiddenSizeL2*hiddenSizeL1), hiddenSizeL2, hiddenSizeL1);
%display_network(W2', 12);                              
%% compute activation from layer 2
[sae2Features] = feedForwardAutoencoder(sae2OptTheta, hiddenSizeL2, ...
                                        hiddenSizeL1, sae1Features);

%% train multi outputs logstic regression                                    
lambda_mr=1e-4;
options_mr.maxIter = 100;
trainLabels=train_labels;
mrModel = mrTrain(hiddenSizeL2, outputSize, lambda_mr, ...
                            sae2Features, trainLabels, options_mr);

saeMultRegOptTheta = mrModel.optTheta(:);

%% fine tuning

% Initialize the stack using the parameters learned
stack = cell(2,1);
inputSize=visibleSize;

stack{1}.w = reshape(sae1OptTheta(1:hiddenSizeL1*inputSize), ...
                     hiddenSizeL1, inputSize);
stack{1}.b = sae1OptTheta(2*hiddenSizeL1*inputSize+1:2*hiddenSizeL1*inputSize+hiddenSizeL1);

stack{2}.w = reshape(sae2OptTheta(1:hiddenSizeL2*hiddenSizeL1), ...
                     hiddenSizeL2, hiddenSizeL1);
stack{2}.b = sae2OptTheta(2*hiddenSizeL2*hiddenSizeL1+1:2*hiddenSizeL2*hiddenSizeL1+hiddenSizeL2);

% Initialize the parameters for the deep model
[stackparams, netconfig] = stack2params(stack);
stackedAETheta = [ saeMultRegOptTheta ; stackparams ];


[stackedAEOptTheta, loss] = minFunc( @(x) stackedAECost(x, ...
      inputSize, hiddenSizeL2, outputSize, netconfig, ...
      lambda, train_input, train_labels), ...
      stackedAETheta, options);
%% save results
%  save ROI_V32_H100_rho_0p01
%% test 
% load test data
% load matFiles/validation_data; 
% [xIsize,yIsize,zIsize]=size(t_I);
% test_input=sampleIMAGES(t_I,patchsize);
% 
% [pyroi] = stackedAEPredict(stackedAEOptTheta, inputSize, hiddenSizeL2, ...
%                           outputSize, netconfig, test_input);
%                                        
% yroi_h=reshape(pyroi,patchsize,patchsize,[]);
% %%
% for k=1:size(yroi_h,3)
%     y1=yroi_h(:,:,k)';
%     y1=imresize(y1,xIsize/patchsize);
%     subplot(4,6,k)
%     imshow(y1)
%     I1=t_I(:,:,k);
%     ymask(:,:,k)=y1.*I1;
% end
% disImgs(ymask,19);
