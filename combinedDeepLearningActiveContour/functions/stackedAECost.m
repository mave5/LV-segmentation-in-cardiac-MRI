function [ cost, grad ] = stackedAECost(theta, inputSize, hiddenSize, ...
                                              outputSize, netconfig, ...
                                              lambda, data, labels)
                                         
% stackedAECost: Takes a trained softmaxTheta and a training data set with labels,
% and returns cost and gradient using a stacked autoencoder model. Used for
% finetuning.
                                         
% theta: trained weights from the autoencoder
% visibleSize: the number of input units
% hiddenSize:  the number of hidden units *at the 2nd layer*
% numClasses:  the number of categories
% netconfig:   the network configuration of the stack
% lambda:      the weight regularization penalty
% data: Our matrix containing the training data as columns.  So, data(:,i) is the i-th training example. 
% labels: A vector containing labels, where labels(i) is the label for the
% i-th training example


%% Unroll softmaxTheta parameter

% We first extract the part which compute the softmax gradient
mrTheta = reshape(theta(1:hiddenSize*outputSize), outputSize, hiddenSize);

% Extract out the "stack"
stack = params2stack(theta(hiddenSize*outputSize+1:end), netconfig);

% You will need to compute the following gradients
mrThetaGrad = zeros(size(mrTheta));
stackgrad = cell(size(stack));
for delta = 1:numel(stack)
    stackgrad{delta}.w = zeros(size(stack{delta}.w));
    stackgrad{delta}.b = zeros(size(stack{delta}.b));
end

cost = 0; % You need to compute this

% You might find these variables useful
M = size(data, 2);
%groundTruth = full(sparse(labels, 1:M, 1));
numCases=M;

%% --------------------------- YOUR CODE HERE -----------------------------
%  Instructions: Compute the cost function and gradient vector for 
%                the stacked autoencoder.
%
%                You are given a stack variable which is a cell-array of
%                the weights and biases for every layer. In particular, you
%                can refer to the weights of Layer d, using stack{d}.w and
%                the biases using stack{d}.b . To get the total number of
%                layers, you can use numel(stack).
%
%                The last layer of the network is connected to the softmax
%                classification layer, softmaxTheta.
%
%                You should compute the gradients for the softmaxTheta,
%                storing that in softmaxThetaGrad. Similarly, you should
%                compute the gradients for each layer in the stack, storing
%                the gradients in stackgrad{d}.w and stackgrad{d}.b
%                Note that the size of the matrices in stackgrad should
%                match exactly that of the size of the matrices in stack.
%


depth = numel(stack);
z = cell(depth+1,1);
a = cell(depth+1, 1);
a{1} = data;

% feed-forward
for layer = (1:depth)
  z{layer+1} = stack{layer}.w * a{layer} + repmat(stack{layer}.b, [1, size(a{layer},2)]);
  a{layer+1} = sigmoid(z{layer+1});
end

% for the output layer we compute cost and gradient 
z_mr = mrTheta * a{depth+1};
hz_mr = sigmoid(z_mr);
temp1 = labels .* log(hz_mr);
temp2 = (1-labels) .* log(1-hz_mr);
cost1 = - sum(sum(temp1+temp2)) ./ numCases;
cost = cost1 + lambda/2 * sum(mrTheta(:) .^ 2);
mrThetaGrad = -1/numCases * (labels - hz_mr) * a{depth+1}' + lambda * mrTheta;

% compute delta3
delta = cell(depth+1);
delta{depth+1} = -(mrTheta' * (labels - hz_mr)) .* a{depth+1} .* (1-a{depth+1});

% compute delta2 and delta1
for layer = (depth:-1:2)
  delta{layer} = (stack{layer}.w' * delta{layer+1}) .* a{layer} .* (1-a{layer});
end

for layer = (depth:-1:1)
  stackgrad{layer}.w = (1/numCases) * delta{layer+1} * a{layer}';
  stackgrad{layer}.b = (1/numCases) * sum(delta{layer+1}, 2);
end


% -------------------------------------------------------------------------

%% Roll gradient vector
grad = [mrThetaGrad(:) ; stack2params(stackgrad)];

end


% You might find this useful
function sigm = sigmoid(x)
    sigm = 1 ./ (1 + exp(-x));
end