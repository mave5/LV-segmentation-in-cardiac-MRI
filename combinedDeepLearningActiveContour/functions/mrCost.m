% softmaxCost.m
function [cost, grad] = mrCost(theta, numOuts, inputSize, lambda, data, labels)
 
% numOuts - the number of number of outputs 
% inputSize - the size N of the input vector
% lambda - weight decay parameter
% data - the N x M input matrix, where each column data(:, i) corresponds to
% a single test set
% labels - an M x 1 matrix containing the labels corresponding for the input data
%
 
% Unroll the parameters from theta
theta = reshape(theta, numOuts, inputSize);
 

numCases = size(data, 2);
 
%groundTruth = full(sparse(labels, 1:numCases, 1));
cost = 0;
 
thetagrad = zeros(numOuts, inputSize);
 
%% ---------- YOUR CODE HERE --------------------------------------
% Instructions: Compute the cost and gradient for softmax regression.
% You need to compute thetagrad and cost.
% The groundTruth matrix might come in handy.
 
[nfeatures, nsamples] = size(data);

zi = theta * data;
%zi=zi+.0005*randn(size(zi));
hzi = 1./(1+exp(-zi));
temp1 = labels .* log(hzi);
temp2 = (1-labels) .* log(1-hzi);

cost = - sum(sum(temp1+temp2)) ./ nsamples;
cost = cost + sum(sum(theta .^ 2)) .* lambda ./ 2;
 
temp3 = labels - hzi;
temp4 = temp3 * data';
thetagrad = - temp4 ./ nsamples;
thetagrad = thetagrad + lambda .* theta;
 
% ------------------------------------------------------------------
% Unroll the gradient matrices into a vector for minFunc
grad = [thetagrad(:)];
end