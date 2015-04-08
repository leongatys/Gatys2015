%Reproduce Figure 2 from Gatys2015
%Figure 2 in the paper was produced with the simulation results stored in
%"binary_signal_fast.mat".
%To generate your own data eg. with different stimuli, use the script found
%in/Figure2/code/ directory
%lgatys 07/04/2015

load('~/Documents/MATLAB/Gatys2015/Figure2/data/binary_signal_fast.mat');
psth = 1000 * psth;
psths = 1000 * psths;
signal = 1000 * signal;
%% Figure 2a and b (b only for binary signals)
fig1 = Figure(1,'size',[100 100]);
range = 2500:2530;
subplot(2,2,[1 2]);plot(range, [psths(range)' psth(range)' signal(range)']);ylim([0 max(signal) + 2]);
legend('PSTH unreliable','PSTH reliable','Presynaptic rate signal')
xlabel('Time in ms')
ylabel('Firing rate in spikes/s')
title('Example signal')
fig1.cleanup
a = min(signal);
b = max(signal);
dbin = .5;
bins = 0:dbin:15;
pc1 = hist(psth(signal==a),bins)/sum((signal==a))/dbin;
pc2 = hist(psth(signal==b),bins)/sum((signal==b))/dbin;
pc1s = hist(psths(signal==a),bins)/sum((signal==a))/dbin;
pc2s = hist(psths(signal==b),bins)/sum((signal==b))/dbin;
subplot(2,2,3);stairs(bins,pc1,'r');hold;stairs(bins,pc2,'r');stairs(bins,pc1s,'b');stairs(bins,pc2s,'b');
xlabel('Firing rate in spikes/s')
ylabel('Probability mass')
legend('Reliable | signal off','Reliable | signal on','Unreliable | signal off','Unreliable | signal on')
title('Conditional firing rate histograms')
fig1.cleanup
%% correlation coefficient
% reliable
corr(signal',psth')
% unreliable 
corr(signal', psths')

%% mutual information (only for binary signals)
dbin = .1;
bins = 0:dbin:15;
% reliable
pc1 = hist(psth(signal==a),bins)/sum((signal==a))/dbin;
pc2 = hist(psth(signal==b),bins)/sum((signal==b))/dbin;
pm = 1/2 * (pc1 + pc2);
margent = -pm .* log2(pm);
margent(isnan(margent)) = 0;
margent = sum(margent);
condent1 = -pc1 .* log2(pc1);
condent1(isnan(condent1)) = 0;
condent1 = sum(condent1);
condent2 = -pc2 .* log2(pc2);
condent2(isnan(condent2)) = 0;
condent2 = sum(condent2);
condent = 1/2 *(condent1+condent2);
mutinf = dbin *  (margent-condent);
% unreliable
pc1s = hist(psths(signal==a),bins)/sum((signal==a))/dbin;
pc2s = hist(psths(signal==b),bins)/sum((signal==b))/dbin;
pms = 1/2 * (pc1s + pc2s);
margents = -pms .* log2(pms);
margents(isnan(margents)) = 0;
margents = sum(margents);
condent1s = -pc1s .* log2(pc1s);
condent1s(isnan(condent1s)) = 0;
condent1s = sum(condent1s);
condent2s = -pc2s .* log2(pc2s);
condent2s(isnan(condent2s)) = 0;
condent2s = sum(condent2s);
condents = 1/2 *(condent1s+condent2s);
mutinfs = dbin * (margents-condents);
[mutinf, mutinfs]

%% mean output firing rate
[mean(psth) mean(psths)]
%% mean and std of membrane potential distribution
[mean(Vt(:)), mean(Vts(:))]
[std(Vt(:)), std(Vts(:))]
%% isi distribution
 isi = []; for i = 1 : params.N_out; isi = [isi diff(find(spTrain(i,:)))];end
 isis = []; for i = 1 : params.N_out; isis = [isis diff(find(spTrains(i,:)))];end
%% coeff. of variation
[std(isi)/mean(isi),std(isis)/mean(isis)]

