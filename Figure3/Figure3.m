% This script reproduces figure 3d from Gatys2015
% lgatys 08/04/15
VarianceRatios = logspace(0,6,500);
result = zeros(1,numel(VarianceRatios));
for i = 1:numel(VarianceRatios)
    result(i) = binary_mutual_information(VarianceRatios(i),1/2);
end
fig = Figure(1, 'size',[50 50]);
semilogx(VarianceRatios, result)
axis square
xlim([0 1000].^2)
set(gca, 'xtick',[1, 10, 100, 1000].^2)
set(gca, 'xticklabel',[1, 10, 100, 1000].^2)
xlabel('Variance ratio')
ylabel('Mutual information (bits)')
fig.cleanup