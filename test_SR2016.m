% Jonathan Garcia
% 2016-12-25

function [evt_log, runtime] = test_SR2016(spks, T)

% Vesicle Priming Dynamics (Hosoi, Sakaba, Neher - 2007)
k_1b = 0.1e-3;  % msec^-1   - rate of Ca-independent vesicle recruitment
k_1C = 1.0e-3;  % 1/(uM*ms) - rate of Ca-dependent vesicle recruitment
Ca_0 = 0.100;   % 100 nM

% Keep track of everything.
el = 0; NS = length(spks);
N_evt_expect = min(1.4*NS + 0.023*T, 3.6*NS + 0.0041*T) + 50;
evt_log = zeros(N_evt_expect, 2);

% Add a spike source.
spkA = SpikingAxon();

% Set up synaptic vesicle pools.
PRM = VesiclePool(7);
RFR = VesiclePool(0);
VAC = VesiclePool(0);

% Create processes to move vesicles between pools.
sync = SynProcess(5.70e-9, 1.62, 12.4, ...
    [5.50e-1, 3.72e-3, 5.32e-7], [0.606, 6.80, 51.3]);
sync.set_facils('pp', 1, 17.6, 11.0, -0.878);
sync.set_facils('Cs', 1, [1.00; 9.38], [2.96; 1.51], [1.68; 7.75]);
sync.set_facils('Ts', 1, [22.5; 45.1], [2.09; 5.97], [-1.51; 0.967]);
sync.set_facils('Cs', 2, [13.1; 103], [135; 292], [1.26; 1.48]);
sync.set_facils('Cs', 3, 213, 80.9, 2.86);

asyn = SynProcess(1.84e-5, 1.61, 4.06, ...
    [1.17e-4, 1.09e-4], [32.8, 130]);
asyn.set_facils('Cs', 1, [25.2; 141], [10.6; 12.4], [0.686; 1.97]);
asyn.set_facils('Cs', 2, 184, 110, 1.26);

refr = SynProcess(1/6.34);

redk = SynProcess(k_1b+k_1C*Ca_0, 1.14, 4.06, ...
    k_1C*[172, 0.295], [0.531, 200]);
redk.set_facils('Ts', 1, 26.0, 19.4, 0.287);
redk.set_facils('Cs', 2, 136, 14.3, 1.34);

% Put everything in easy-to-access arrays.
processes = [sync; asyn; refr; redk];
responses = {@S_rel, @A_rel, @end_rfr, @redock};
next_T = Inf+zeros(4, 1); time = 0;

% Ensure that this keeps updated on all
% the next event times for each process.
for j = 1:length(processes)
    addlistener(processes(j), 'EventUpdate', ...
        @(src, event) update_next(j, event));
end

% Have all processes listen for spiking events, and connect
% all the pools to the processes that rely on them.
addlistener(spkA, 'Spike', ...
    @(src, event) sync.respond2spike(src, PRM, event));
addlistener(PRM, 'PoolUpdate', ...
    @(src, event) sync.respond2pool(spkA, src, event));
sync.new_evt_time(PRM.pool_size, time, -Inf);

addlistener(spkA, 'Spike', ...
    @(src, event) asyn.respond2spike(src, PRM, event));
addlistener(PRM, 'PoolUpdate', ...
    @(src, event) asyn.respond2pool(spkA, src, event));
asyn.new_evt_time(PRM.pool_size, time, -Inf);

addlistener(spkA, 'Spike', ...
    @(src, event) refr.respond2spike(src, RFR, event));
addlistener(RFR, 'PoolUpdate', ...
    @(src, event) refr.respond2pool(spkA, src, event));
refr.new_evt_time(RFR.pool_size, time, -Inf);

addlistener(spkA, 'Spike', ...
    @(src, event) redk.respond2spike(src, VAC, event));
addlistener(VAC, 'PoolUpdate', ...
    @(src, event) redk.respond2pool(spkA, src, event));
redk.new_evt_time(VAC.pool_size, time, -Inf);

% Run simulation.
tic;
ns = 1;     % index of current spike
while time < T
    [time, evt] = min(next_T);
    if time < new_spike(ns)
        responses{evt}(time);
%         fprintf('(PRM: %d, RFR: %d, VAC: %d)\n', ...
%             PRM.pool_size, RFR.pool_size, VAC.pool_size);
        el = el + 1;
        evt_log(el,:) = [time, evt];
    else
        time = new_spike(ns);
        spkA.spike(time);
        ns = ns + 1;
%         fprintf('Spike: %.4f ms\n', time);
        el = el + 1;
        evt_log(el,:) = [time, 0];
    end
end
runtime = toc;
evt_log = evt_log(1:el,:);

% Helper functions for handling next-event-time updates.

    function update_next(ID, event)
        next_T(ID) = event.new_time;
    end

    function nt = new_spike(ns)
        if ns > NS
            nt = Inf;
        else
            nt = spks(ns);
        end
    end

% Define what happens during events for each process.

    function S_rel(new_time)
        VAC.update(+1, sync, new_time);
        RFR.update(+PRM.pool_size-1, sync, new_time);
        PRM.update(-PRM.pool_size, sync, new_time);
%         fprintf('S_rel: %.4f ms\n', new_time);
    end

    function A_rel(new_time)
        VAC.update(+1, asyn, new_time);
        RFR.update(+PRM.pool_size-1, asyn, new_time);
        PRM.update(-PRM.pool_size, asyn, new_time);
%         fprintf('A_rel: %.4f ms\n', new_time);
    end

    function end_rfr(new_time)
        RFR.update(-1, refr, new_time);
        PRM.update(+1, refr, new_time);
%         fprintf('end_rfr: %.4f ms\n', new_time);
    end

    function redock(new_time)
        VAC.update(-1, redk, new_time);
        PRM.update(+1, redk, new_time);
%         fprintf('redock: %.4f ms\n', new_time);
    end

end