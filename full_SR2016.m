% Jonathan Garcia
% 2017-04-26

function evt_log = full_SR2016(spks, T, nAZ, ...
    nRRP, nREC, sr0, S0s, SFs, ar0, A0s, AFs, ...
    facil_redock, refr_TF, RID_TF, FDR_TF, deplete_TF)

% Set default flags.
if isempty(facil_redock)
    facil_redock = false;   % constant recycling rate
end
if isempty(refr_TF)
    refr_TF = true;         % yes to refractory period
end
if isempty(RID_TF)
    RID_TF = false;         % no release-independent depression
end
if isempty(FDR_TF)
    FDR_TF = false;         % no frequency-dependent recovery
end
if isempty(deplete_TF)
    deplete_TF = true;      % yes to depletion upon release
end

% Keep track of everything.
el = 0; NS = length(spks);
% N_evt_expect = ceil(min(1.4*NS + 0.023*T, 3.6*NS + 0.0041*T) + 50);
evt_log = zeros(10 * NS, 3);

% Add a spike source.
spkA = SpikingAxon();

% Set up synaptic vesicle pools.
PRM(nAZ,1) = VesiclePool;
RFR(nAZ,1) = VesiclePool;
DEP(nAZ,1) = VesiclePool;
REC = VesiclePool(nREC);
for az = 1:nAZ
    PRM(az) = VesiclePool(nRRP(az), nRRP(az));
    RFR(az) = VesiclePool(0);
    DEP(az) = VesiclePool(0);
end

% Create processes to move vesicles between pools.
sync(nAZ,1) = SynProcess;
asyn(nAZ,1) = SynProcess;
refr(nAZ,1) = SynProcess;
kRID(nAZ,1) = SynProcess;
kFDR(nAZ,1) = SynProcess;
redk(nAZ,1) = SynProcess;
for az = 1:nAZ
    % synchronous release at all active zones
    sync(az) = SynProcess(sr0, S0s{az});
    sync(az).set_facils(1, SFs{1,1}, SFs{1,2}, SFs{1,3});
    sync(az).set_facils(2, SFs{2,1}, SFs{2,2}, SFs{2,3});
    sync(az).set_facils(3, SFs{3,1}, SFs{3,2}, SFs{3,3});
%     sync(az).set_facils(4, SFs{4,1}, SFs{4,2}, SFs{4,3});
    
    % asynchronous release at all active zones
    asyn(az) = SynProcess(ar0, A0s{az});
    asyn(az).set_facils(1, AFs{1,1}, AFs{1,2}, AFs{1,3});
    asyn(az).set_facils(2, AFs{2,1}, AFs{2,2}, AFs{2,3});
%     asyn(az).set_facils(3, AFs{3,1}, AFs{3,2}, AFs{3,3});
    
    % end of release refractory period at each active zone
    refr(az) = SynProcess(1/6.34);
    
    % release-independent depression and frequency-dependent recovery
    if RID_TF
        kRID(az) = SynProcess(0.00, ...
            {{0.400, 10.0, 0.500, 4.00, 0.500}});
        kRID(az).set_facils(1, 10.0, 10.0, 1.00);
    end
    
    if FDR_TF
        kFDR(az) = SynProcess(1/2000, ...
            {{2.40e-1, 320., 0.500, 4.00, 0.500}});
        kFDR(az).set_facils(1, 320., 10.0, 2.00);
    else
        kFDR(az) = SynProcess(1/2000);  % constant recovery
    end
    
    % vesicle redock rate from recycling pool to each active zone
    if facil_redock
        % Vesicle Priming Dynamics (Hosoi, Sakaba, Neher - 2007)
        k_1b = 0.1e-3;  % msec^-1   - rate of Ca-independent recruitment
        k_1C = 1.0e-3;  % 1/(uM*ms) - rate of Ca-dependent recruitment
        Ca_0 = 0.100;   % 100 nM
        
        redk(az) = SynProcess(k_1b+k_1C*Ca_0, ...
            {{5.00e-2, 1.00e-1, 1.75e+0, 3.41e+0, 1.75e-1}, ...
             {6.00e-2, 8.00e+1, 5.25e-1, 1.00e+1, 4.50e+0}, ...
             {5.00e-2, 1.00e+3, 1.50e-1, 5.00e+1, 1.00e+1}});
        redk(az).set_facils(1, 1.00e+1, 1.00e+1, 1.00e+0);
        redk(az).set_facils(2, 2.00e+2, 1.25e+1, 2.50e+0);
    else
        redk(az) = SynProcess(1/2800);
    end
end

% Put everything in easy-to-access arrays.
processes = [sync; asyn; refr; kRID; kFDR; redk];
responses = {@S_rel, @A_rel, @end_rfr, @depress, @recover, @redock};
next_T = Inf+zeros(length(processes), 1);
time = 0;

% Ensure that this keeps updated on all
% the next event times for each process.
for j = 1:length(processes)
    addlistener(processes(j), 'EventUpdate', ...
        @(src, event) update_next(j, event));
end

% Have all processes listen for spiking events, and connect
% all the vesicle pools to the processes that rely on them.
for az = 1:nAZ
    connect(sync(az), spkA, PRM(az), REC);
    connect(asyn(az), spkA, PRM(az), REC);
    connect(refr(az), spkA, RFR(az), PRM(az));
    connect(kRID(az), spkA, PRM(az), DEP(az));
    connect(kFDR(az), spkA, DEP(az), PRM(az));
    connect(redk(az), spkA, REC, PRM(az));
end

% Run simulation.
% tic;
ns = 1;     % index of current spike
while time < T
    [time, evt] = min(next_T);
    if time < new_spike(ns)
        [id, az] = proc_ID(evt);
        responses{id}(time, az);
%         fprintf('(PRM: %d, RFR: %d, VAC: %d)\n', ...
%             PRM.pool_size, RFR.pool_size, VAC.pool_size);
        el = el + 1;
        evt_log(el,:) = [time, id, az];
    else
        time = new_spike(ns);
        spkA.spike(time);
        ns = ns + 1;
%         fprintf('Spike: %.4f ms\n', time);
        el = el + 1;
        evt_log(el,:) = [time, 0, 0];
    end
end
% runtime = toc;
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

    function S_rel(new_time, az)
        if deplete_TF
            transfer(sync(az), 1, PRM(az), REC, new_time, false);
        else
            transfer(sync(az), 0, PRM(az), REC, new_time, false);
        end
        if refr_TF
            transfer(sync(az), PRM(az).pool_size, ...
                PRM(az), RFR(az), new_time, true);
        end
%         fprintf('S_rel: %.4f ms\n', new_time);
    end

    function A_rel(new_time, az)
        if deplete_TF
            transfer(asyn(az), 1, PRM(az), REC, new_time, false);
        else
            transfer(asyn(az), 0, PRM(az), REC, new_time, false);
        end
        if refr_TF
            transfer(asyn(az), PRM(az).pool_size, ...
                PRM(az), RFR(az), new_time, true);
        end
%         fprintf('A_rel: %.4f ms\n', new_time);
    end

    function end_rfr(new_time, az)
        transfer(refr(az), 1, RFR(az), PRM(az), new_time, true);
%         fprintf('end_rfr: %.4f ms\n', new_time);
    end

    function depress(new_time, az)
        transfer(kRID(az), 1, PRM(az), DEP(az), new_time, true);
    end

    function recover(new_time, az)
        transfer(kFDR(az), 1, DEP(az), PRM(az), new_time, true);
    end

    function redock(new_time, az)
        transfer(redk(az), 1, REC, PRM(az), new_time, false);
%         fprintf('redock: %.4f ms\n', new_time);
    end

% Build a connection of a process from one vesicle pool to another.
    function connect(process, spk_src, ves_src, ves_dst)
        % Listen for the spikes.
        addlistener(spk_src, 'Spike', @(~, event) ...
            process.respond2spk(spk_src, ves_src, ves_dst, event));
        
        % Listen for changes to the source pool.
        addlistener(ves_src, 'PoolUpdate', @(~, event) ...
            process.respond2src(spk_src, ves_src, ves_dst, event));
        
        % Listen for changes to the destination pool.
        addlistener(ves_dst, 'PoolUpdate', @(~, event) ...
            process.respond2dst(spk_src, ves_src, ves_dst, event));
        
        % Initialize first event time.
        N_ves = min(ves_src.pool_size, ves_dst.vacancies);
        process.new_evt_time(N_ves, 0, [-Inf, -Inf]);
    end

% Transfer vesicles from one pool/state to another.
    function transfer(process, N_ves, ves_src, ves_dst, new_time, share_TF)
        ves_src.update(-N_ves, process, new_time, share_TF);
        ves_dst.update(+N_ves, process, new_time, share_TF);
    end

% Convert (non-spiking) event index into an event-type/active-zone pair.
    function [id, az] = proc_ID(evt)
        id = floor((evt - 1) / nAZ) + 1;
        az = mod(evt - 1, nAZ) + 1;
    end

end

% for az = 1:nAZ
%     addlistener(spkA, 'Spike', ...
%         @(~, event) sync(az).respond2spk(spkA, PRM(az), REC, event));
%     addlistener(PRM(az), 'PoolUpdate', ...
%         @(~, event) sync(az).respond2src(spkA, PRM(az), REC, event));
%     addlistener(REC, 'PoolUpdate', ...
%         @(~, event) sync(az).respond2dst(spkA, PRM(az), REC, event));
%     sync(az).new_evt_time(PRM(az).pool_size, 0, [-Inf, -Inf]);
% 
%     addlistener(spkA, 'Spike', ...
%         @(src, event) asyn(az).respond2spk(spkA, PRM(az), REC, event));
%     addlistener(PRM(az), 'PoolUpdate', ...
%         @(src, event) asyn(az).respond2src(spkA, PRM(az), REC, event));
%     addlistener(REC, 'PoolUpdate', ...
%         @(src, event) asyn(az).respond2dst(spkA, PRM(az), REC, event));
%     asyn(az).new_evt_time(PRM(az).pool_size, 0, [-Inf, -Inf]);
% 
%     addlistener(spkA, 'Spike', ...
%         @(src, event) refr(az).respond2spk(spkA, RFR(az), PRM(az), event));
%     addlistener(RFR(az), 'PoolUpdate', ...
%         @(src, event) refr(az).respond2src(spkA, RFR(az), PRM(az), event));
%     addlistener(PRM(az), 'PoolUpdate', ...
%         @(src, event) refr(az).respond2dst(spkA, RFR(az), PRM(az), event));
%     refr(az).new_evt_time(RFR(az).pool_size, 0, [-Inf, -Inf]);
% 
%     addlistener(spkA, 'Spike', ...
%         @(src, event) redk(az).respond2spk(spkA, REC, PRM(az), event));
%     addlistener(REC, 'PoolUpdate', ...
%         @(src, event) redk(az).respond2src(spkA, REC, PRM(az), event));
%     addlistener(PRM(az), 'PoolUpdate', ...
%         @(src, event) redk(az).respond2dst(spkA, REC, PRM(az), event));
%     redk(az).new_evt_time(REC.pool_size, 0, [-Inf, -Inf]);
% end