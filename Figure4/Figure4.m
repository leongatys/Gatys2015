%This script reproduces Figure 4 in Gatys2015
%lgatys 08/04/15
%% load data for reliable synapses and bring it in the right form
load('~/Documents/MATLAB/Gatys2015/Figure4/data/gradient_search_reliable.mat')
nSample = 2 * round(logspace(1, log10(25000), 1000) / 2);
nSample = 2 * nSample;
data_reliable = struct2cell(data_reliable);
params_reliable = squeeze(cell2mat(data_reliable(1, 1, :)))';
mi_reliable = -squeeze(cell2mat(data_reliable(2, 1, :)));
mi_reliable = flipud(mi_reliable);

%% load data for unreliable synapses and bring it in the right form
load('~/Documents/MATLAB/Gatys2015/Figure4/data/gradient_search_unreliable.mat')
nSample = 2 * round(logspace(1, log10(25000), 1000) / 2);
nSample = 2 * nSample;
data_unreliable = struct2cell(data_unreliable);
params_unreliable = squeeze(cell2mat(data_unreliable(1, 1, :)))';
mi_unreliable = -squeeze(cell2mat(data_unreliable(2, 1, :)));
mi_unreliable = flipud(mi_unreliable);
%positions and weights of delta peaks of input signal
pos_unreliable =  flipud([1e-6  * ones(length(params_unreliable), 1) params_unreliable(:, 1 : 3) .1 * ones(length(params_unreliable), 1)]);
wei_unreliable = flipud([params_unreliable(:, 4 : 7) 1 - sum(params_unreliable(:, 4 : 7),2)]);
%merge peaks that are on top of each other into one for plotting
ndx1 = abs(pos_unreliable(:, 3)./pos_unreliable(:, 2) - 1) < 0.1;
ndx2 = abs(pos_unreliable(:, 3)./pos_unreliable(:, 4) - 1) < 0.1;
ndx3 = abs(pos_unreliable(:, 4)./pos_unreliable(:, 2) - 1) < 0.1;
wei_unreliable(ndx1, 2) = wei_unreliable(ndx1, 2) + wei_unreliable(ndx1, 3);
wei_unreliable(ndx3, 2) = wei_unreliable(ndx3, 2) + wei_unreliable(ndx3, 4);
wei_unreliable(ndx2, 4) = wei_unreliable(ndx2, 3) + wei_unreliable(ndx2, 4);
wei_unreliable(ndx1, 3) = NaN;
wei_unreliable(ndx2, 3) = NaN;
wei_unreliable(ndx3, 4) = NaN;
pos_unreliable(isnan(wei_unreliable)) = NaN;
ndx = find(wei_unreliable == 0);
pos_unreliable(ndx) = NaN;
wei_unreliable(ndx) = NaN;
%% Plot optimal signal values and weights and corresponding channel capacity
fig = Figure(1, 'size',[100 200]);

subplot(3,1,1)
loglog(nSample, pos_unreliable);
xlim([nSample(1) nSample(end)])
ylim([1e-7 1])
set(gca, 'xtick',[100 1000 10000 50000])
set(gca,'xticklabel',[100 1000 10000 50000])
set(gca, 'ytick',[1e-6 1e-4 1e-3 1e-2 1e-1])
set(gca,'yticklabel',[0 0.1 1 10 100])
ylabel('Firing rate (spikes/s)')

subplot(3,1,2)
semilogx(nSample, wei_unreliable);xlim([0 50000])
ylim([0,.6])
xlim([nSample(1) nSample(end)])
set(gca, 'xtick',[100 1000 10000 50000])
set(gca,'xticklabel',[100 1000 10000 50000])
ylabel('Probability mass')

subplot(3,1,3)
semilogx(nSample, mi_unreliable,'k');hold
semilogx(nSample, mi_reliable,'color',[1 1 1]/2)
xlim([nSample(1) nSample(end)])
ylim([0 1.1])
set(gca, 'xtick',[100 1000 10000 50000])
set(gca,'xticklabel',[100 1000 10000 50000])
ylabel('Channel Capacity (bits)')
xlabel('Population size (neurons)')
fig.cleanup

%% crossections at 5,000 and 50,000 neurons
pos5k = pos_unreliable(end-294,:);
pos5k = pos5k(~isnan(pos5k));
wei5k = wei_unreliable(end-294,:);
wei5k = wei5k(~isnan(wei5k));

pos50k = pos_unreliable(end,:);
pos50k = pos50k(~isnan(pos50k));
wei50k = wei_unreliable(end,:);
wei50k = wei50k(~isnan(wei50k));

fig3 = Figure('size',[30 50]);
subplot(2,1,1)
bar(log(pos5k),wei5k,.5)
set(gca,'xtick',log(pos5k))
set(gca,'xticklabel',round(10000*pos5k)/10)
set(gca,'xlim',[log(min(pos5k))-2 log(max(pos5k))+2])
set(gca,'ylim',[0 0.5])
set(gca,'ytick',(0:.1:.5))
axis square

subplot(2,1,2)
bar(log(pos50k),wei50k,.75)
set(gca,'xtick',log(pos50k))
set(gca,'xticklabel',round(10000*pos50k)/10)
set(gca,'xlim',[log(min(pos50k))-2 log(max(pos50k))+2])
set(gca,'ylim',[0 0.5])
set(gca,'ytick',(0:.1:.5))
xlabel('Firing rate (spikes/s)')
ylabel('Probability mass')
axis square

fig3.cleanup