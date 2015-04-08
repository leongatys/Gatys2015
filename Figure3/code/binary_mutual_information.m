function[mutualInformation] = binary_mutual_information(ratio,weight)

% Computes the mutual information between a binary signal with equal
% probability mass on both signal states.
% numerically approximates the conditional postsynaptic current distribution by a
% normal distribution
% ratio is the ratio of the conditional variances
% weight is the probability mass on the smaller signal state
% lgatys 07/04/2015

weightSignal  = [weight, 1-weight]; %the probability mass on the signal values
margEntropy = 0;
condEntropy = 0;
% Compute conditional response distributions J|R=r
d = 1; %how finely to sample the normal distribution
range = 10000; %range of sampling the normal distribution
condResponseDist = zeros(range/d + 1, 2);
condResponseDist(:,1) = d * normpdf(0:d:range,0,1);
condResponseDist(:,2) = d * normpdf(0:d:range,0,sqrt(ratio));

%Compute the marginal entropy of resulting process J
pJ = condResponseDist * weightSignal';
mEntVec = - pJ .* log2(pJ);
mEntVec(isnan(mEntVec)) = 0;
margEntropy = 2 * sum(mEntVec) - mEntVec(1);

% Computing the conditional entropy of J given R
for r = 1 : 2
    J = 1;
    pCondJ = condResponseDist(J, r);
    dcondEntropy = -pCondJ * log2(pCondJ);
    pCondJ = condResponseDist(2 : end, r);
    ndx = pCondJ > 0;
    dcondEntropy = dcondEntropy - 2 * pCondJ(ndx)' * log2(pCondJ(ndx));
    condEntropy = condEntropy + weightSignal(r) * dcondEntropy;
end
mutualInformation = margEntropy - condEntropy;
end