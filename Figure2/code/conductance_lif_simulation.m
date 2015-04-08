% This script implements the simulation of a conductance-based leaky
% integrate and fire neuron as found in Gatys2015.
% The neuron receives balanced presynaptic input and synaptic transmission
% is either reliable or unreliable.
% lgatys 07/04/15
%% presynaptic rate signal
T = 1000; params.T = T;%total simulation time in ms
switchTime = 1;%time in ms after which the signal is redrawn
maxRate = .01;%maximum input firing rate in spikes/ms
signal = maxRate * (kron(rand(1,T),ones(switchTime,1))<1/2);%binary signal
% signal = maxRate * kron(rand(1,T),ones(switchTime,1));%uniform signal
signal = signal(1:T);

%% simulation
%parameter
dt = 1; params.dt = dt; % timestep in ms
tau_p = 22; params.tau_leak = tau_p; %leak time constant in ms
Vth = -55.3; params.Vthreshold = Vth; %spiking threshold in mV for reliable synapses
Vths = -52.95; params.Vthreshold_s =Vths; %spiking threshold in mV for unreliable synapses
Vr = -70; params.Vreset = Vr; %reset voltage in mV
V0 = -80; params.Vrest = V0;% resting potential in mV
Ve = 0; params.Vreversalex = Ve;%ex.reversal potential in mV
Vi = -75; params.Vreversalin = Vi;%in. reversal potenatial in mV
ge = 1.* .0032; params.ge = ge;%scaling factors of the impact of a synaptic event on the postsynaptic conductances 
gi = 1.* .064; params.gi = gi;
%c =0.35; params.membraneCapacitance = c;% membrane capacitance in nF not
%needed in simulation.
N_ex = 4000; params.N_ex = N_ex;%ex input population for each neuron
N_in = 1000; params.N_in = N_in;%inhib input population for each neuron
N_out = 10000; params.N_out = N_out;%total output population to estimate psth

spTrain = zeros(N_out,T);%output spike train reliable synapses
spTrains = zeros(N_out,T);%output spike train unreliable synapses
Vt = zeros(N_out,T);%postsynaptic membrane potential trace reliable synapses
Vts = zeros(N_out,T);%postsynaptic membrane potential trace unreliable synapses
parfor n = 1 : N_out;
    tic
    V = -65; % starting membrane voltage
    Vs = -65;
    % compute input spike trains with reliable and unreliable synapses
    exSp = zeros(1,T);
    inSp = zeros(1,T);
    exSpS = zeros(1,T);
    inSpS = zeros(1,T);
    for k = 1 : N_ex
        ex = rand(1, T) < signal;
        exSp = exSp + ex;
        exSpS(ex) = exSpS(ex) + poissrnd(1,1,sum(ex));
    end
    for k = 1 : N_in
        in = rand(1, T) < signal;
        inSp = inSp + in;
        inSpS(in) = inSpS(in) + poissrnd(1,1,sum(in));
    end
    
    
    % Simulate conductance based Lif Neuron
    % diff equation: 
    % dV/dt = -(V-V_p)/tau_p - (g_e (V-V_e) + g_i (V-V_i)) + epsilon
    % ref ensures the refractory period of 2ms
    ref = 0;
    refs = 0;
    for i = 1 : T
        
        eps = 2 * randn(1);
        if ref <= 0
            V = V  - (V-V0)/tau_p - (ge * (V - Ve) * exSp(i) + gi * (V - Vi) * inSp(i)) + eps;
        end
        
        if refs <= 0
            Vs = Vs  - (Vs-V0)/tau_p - (ge * (Vs - Ve) * exSpS(i) + gi * (Vs - Vi) * inSpS(i)) + eps;
        end
        
        ref = ref - 1;
        refs = refs - 1;
        
        if V >= Vth
            spTrain(n,i) = 1;
            V = Vr;
            ref = 2;
        end
        if Vs >= Vths
            spTrains(n,i) = 1;
            Vs = Vr;
            refs = 2;
        end
        
        
        Vt(n,i) = V;
        Vts(n,i) = Vs;
        
    end
    toc
end
psth = mean(spTrain);
psths = mean(spTrains);

%save('MYDIR/Results.mat','params','psth','psths','Vt','Vts','spTrain','spTrains','signal')
