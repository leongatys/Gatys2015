%This script produces data with unreliable synapses for Figure 4 in Gatys2015
%For different number of input neurons it maximizes the mutual information between a signal consisting of 5 delta
%peaks at positions between 0 and 0.1 and the postsynaptic current. 
%The maximization is performed iteratively - first the maximum mutual
%information is found for the largest value, then this result is used as a
%starting point for the gradient descent at the second largest value etc.
%lgatys 08/04/15
nSample = 2 * round(logspace(1, log10(25000), 1000) / 2); %number of neurons at which to find the maximum 
nSample = fliplr(nSample);
startparams = [0.0007    0.0001    0.0034    0.3692    0.0799    0.1619    0.1102]; %these params where found with gridsearch and gradient descent for 50000 input neurons but any reasonable init. will do
for i = 1 : length(nSample);
    i
    tic
    numNeurons = nSample(i);
    mutInformation_gd = @(arg)-mutual_information_unreliable(numNeurons, arg);
    options = optimset('TolX', 1e-10, 'TolFun', 1e-10, 'MaxFunEvals', 7000,'UseParallel','always');
    if i == 1
        [data_unreliable(i).par data_unreliable(i).mi] = fmincon(mutInformation_gd, startparams, [0 0 0 1 1 1 1], 1, [], [], [0 0 0 0 0 0 0], [.1 .1 .1 1 1 1 1], [], options);
    else
        [data_unreliable(i).par data_unreliable(i).mi] = fmincon(mutInformation_gd, data_unreliable(i - 1).par, [0 0 0 1 1 1 1], 1, [], [], [0 0 0 0 0 0 0], [.1 .1 .1 1 1 1 1], [], options);
        
    end
    
    toc
end
% save('~/Documents/MATLAB/Gatys2015/Figure4/data/gradient_search_unreliable.mat','data_unreliable')