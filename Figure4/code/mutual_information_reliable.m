function[mutualInformation, margEntropy] = mutual_information_reliable(numNeurons, arg)

% Computes the mutual information between signal that is a collection of
% delta peaks at positions posSignal with weighting weightSignal and resulting 
% process J in a model with 2 * numNeurons neurons.
% For N delta peaks the first N - 2 entries of arg are the positions of
% the N - 2 peaks in between the highest and lowest peak. The last N - 1
% entries of arg are the weights on the lowest N - 1 peaks and the weight
% on the highest peak follows from that (total probability = 1).
% Adds standard normal noise on the PSC
% lgatys 08-04-2015

highPos = 0.1; %two delta peaks are fixed on highest and lowest possible signal value
lowPos = 0;
posSignal = [lowPos arg(1 : floor(numel(arg) / 2)) highPos];
weightSignal = arg(floor(numel(arg) / 2) + 1 : end);
weightSignal  = [weightSignal (1 - sum(weightSignal))];
condEntropy = 0;
% Compute conditional response distributions J|R=r
binoDist = zeros(numNeurons + 1, numel(posSignal));
condResponseDist = zeros(2 * numNeurons + 1, numel(posSignal));
noiseDist = normpdf(-10 : 10, 0, 1);
for i = 1 : numel(posSignal)
    binoDist(:, i) = binopdf((0:numNeurons), numNeurons, abs(posSignal(i)));
end
for i = 1 : numel(posSignal)
    condResponseDist(:,i) = conv(binoDist(:,i), flipud(binoDist(:,i)));
    %Add GWN on the PSC
    condResponseDist(:,i) = conv(condResponseDist(:, i), noiseDist,'same');
end
condResponseDist = condResponseDist(numNeurons + 1 : end, :); %use symmetry to save computation time

%Compute the marginal entropy of resulting process J
pJ = condResponseDist * weightSignal';
mEntVec = - pJ .* log2(pJ);
mEntVec(isnan(mEntVec)) = 0;
margEntropy = 2 * sum(mEntVec) - mEntVec(1);

% Computing the conditional entropy of J given R
for r = 1 : length(posSignal)
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