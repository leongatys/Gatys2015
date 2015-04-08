function[mutualInformation, margEntropy] = mutual_information_unreliable(numNeurons, arg)

% Computes the mutual information between signal that is a collection of
% delta peaks at positions posSignal with weighting weightSignal and
% resulting process J in a model with 2 * numNeurons neurons.
% Synaptic transmission is noisy and modelled as a poisson process with
% mean 1.
% For N delta peaks the first N - 2 argument entries are the positions of
% the N - 2 peaks in between the highest and lowest peak. The last N - 1
% argument entries are the weights on the lowest N - 1 peaks and the weight
% on the highest peak follows from that (total probability = 1 ).
% Adds standard normal noise on the PSC
% lgatys 08-04-2015
highPos = 0.1; %two delta peaks are fixed on highest and lowest possible signal value
lowPos = 0;
posSignal = [lowPos arg(1 : floor(numel(arg) / 2)) highPos];
weightSignal = arg(floor(numel(arg) / 2) + 1 : end);
weightSignal  = [weightSignal (1 - sum(weightSignal))];
condEntropy = 0;
% Compute conditional response distributions J|R=r
neuronDist = zeros(11, numel(posSignal));
condResponseDist = zeros(2 * 10 * numNeurons + 1, numel(posSignal));
noiseDist = normpdf(-10 : 10, 0, 1); 
for i = 2 : numel(posSignal)
    %postsynaptic current distribution from a single neuron
    neuronDist(:, i) = posSignal(i) .* poisspdf(0:10,1);
    neuronDist(1, i) = neuronDist(1, i) + (1 - posSignal(i));
    %postsynaptic current distribution from whole presynaptic population
    popDist = zeros(1, 10 * numNeurons + 1);
    popDist(1, 1 : 10 + 1) = neuronDist(:, i);
    popDist = abs(ifft(fft(popDist).^(numNeurons - 1)))';
    condResponseDist(:, i) = conv(popDist, flipud(popDist));
    %Add GWN on the PSC
    condResponseDist(:,i) = conv(condResponseDist(:, i), noiseDist,'same');
end
condResponseDist = condResponseDist(10 * numNeurons + 1 : end, :);
condResponseDist(1 : ceil(length(noiseDist)/2), 1) = noiseDist(ceil(length(noiseDist/2)/2 : end));%for r=0 PSC is only the additive noise
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