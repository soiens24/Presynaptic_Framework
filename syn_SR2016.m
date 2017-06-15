% Jonathan Garcia
% 2016-04-21

function [rels, Prel, p, history, vp_init] = syn_SR2016(spks, p, nAZ, T)
% Presynaptic vesicle release model:
% Event-driven vesicle release based on spike times and vesicle pool sizes.

% Set default parameters.
if nargin < 4
    T = [];
if nargin < 3
    nAZ = [];
if nargin < 2
    p = [];
end
end
end
if isempty(p)
    p = struct;
    
    % default vesicle pool sizes
    p.RRP0 = ...    % RPP <-> Readily Releasable Pool -> KAR/EMT
        [0; ...     % vacant spots
         0; ...     % inactivated/depressed vesicles
         0; ...     % in post-release refractory period
         7];        % ready to release
    p.KAR0 =  0;    % RRP -> Kiss-And-Run pool -> PMP
    p.EMT0 =  0;    % RRP -> EMpTy pool -> PMP
    p.PMP0 =  0;    % EMT -> PreMature Priming pool -> RPP
    p.RPP0 = 10;    % PMP -> Readily Priming Pool <-> RRP
    
    % synchronous release rate profile parameters
    p.sr0 = 5.70e-9;    % ms^-1 - steady-state release rate
    p.S_C = ...         % ms^-1 - triple exponential coefficients
        [2.83e-1; ...       % fast
         2.10e-3; ...       % medium
         2.30e-7];          % slow
    p.S_T = ...         % msec  - triple exponential time constants
        [7.15e-1; ...       % fast
         6.70e+0; ...       % medium
         5.02e+1];          % slow
    p.S_h = 1.58e+0;    % msec  - initial ramp 0.5-rise time
    p.S_p = 1.41e+1;    % power on sigmoid
%     p.sr0 = 5.7006e-9;  % ms^-1 - steady-state release rate
%     p.S_C = ...         % ms^-1 - triple exponential coefficients
%         [2.2726e-1; ...     % fast
%          2.9233e-3; ...     % medium
%          8.8341e-7];        % slow
%     p.S_T = ...         % msec  - triple exponential time constants
%         [7.5347e-1; ...     % fast
%          6.5877e+0; ...     % medium
%          5.0030e+1];        % slow
%     p.S_h = 1.5827e+0;  % msec  - initial ramp 0.5-rise time
%     p.S_p = 1.3870e+1;  % power on sigmoid
    
    % synchronous release rate facilitation parameters
    p.S_F_c1 = [ 3.83e-01,  1.83e-03; ...	
                 5.88e+00,  1.85e+02; ...
                 8.83e+01,  9.05e+01; ...
                +1,        +1];
    p.S_F_t1 = [ 8.78e-01,  8.96e-07; ...
                 1.92e+01,  6.80e+01; ...
                 1.52e+00,  4.34e+05; ...
                -1,        +1];
    p.S_F_c2 = [ 5.37e-22,  1.55e-05; ...
                 9.11e+00,  1.50e+02; ...
                 1.05e+22,  1.60e+05; ...
                +1,      +1];
    p.S_F_c3 = [ 1.99e-07; ...
                 1.31e+02; ...
                 1.89e+08; ...
                +1];
    p.S_F_p0 = [ 1.20e-01; ...
                 1.52e+01; ...
                 9.29e+00; ...
                -1];
%     p.S_F_c1 = [ 2.4853e+1,  1.8429e+0; ...     % ha, hb
%                  5.9536e+0,  1.0693e+2; ...     % Ta, Tb
%                  9.3303e-1,  4.9916e-1];        % pa, pb
%     p.S_F_t1 = [ 7.7334e+0,  1.0346e+0; ...     % ha, hb
%                  1.4570e+1,  3.2883e+1; ...     % Ta, Tb
%                 -5.1323e-1,  9.9736e-1];        % pa, pb
%     p.S_F_c2 = [ 4.8522e-1,  5.0702e-1; ...     % ha, hb
%                  1.1729e+2,  1.4219e+1; ...     % Ta, Tb
%                  2.6154e+0,  2.4995e+0];        % pa, pb
%     p.S_F_c3 = [ 3.9976e-1; ...                 % h
%                  2.2583e+2; ...                 % T
%                  4.4338e+0];                    % p
%     p.S_F_p0 = [ 4.2703e+0; ...                 % h
%                  1.0909e+1; ...                 % T
%                 -5.0389e-1];                    % p
    
    % asynchronous release rate profile parameters
    p.ar0 = 1.84e-5;    % ms^-1 - steady-state release rate
    p.A_C = ...         % ms^-1 - double exponential coefficients
        [1.16e-4; ...       % fast
         8.42e-5];          % slow
    p.A_T = ...         % msec  - double exponential time constants
        [3.54e+1; ...       % fast
         1.46e+2];          % slow
    p.A_h = 1.57e+0;    % msec  - initial ramp half-rise time
    p.A_p = 2.96e+0;    % power on sigmoid
%     p.ar0 = 1.8357e-5;  % ms^-1 - steady-state release rate
%     p.A_C = ...         % ms^-1 - double exponential coefficients
%         [1.5486e-4; ...     % fast
%          9.8871e-5];        % slow
%     p.A_T = ...         % msec  - double exponential time constants
%         [3.3806e+1; ...     % fast
%          1.2898e+2];        % slow
%     p.A_h = 1.5820e+0;  % msec  - initial ramp half-rise time
%     p.A_p = 2.8770e+0;  % power on sigmoid
    
    % asynchronous release rate facilitation parameters
    p.A_F_c1 = [ 1.41e-07,  9.08e-08; ...
                 2.40e+02,  4.39e+01; ...
                 1.23e+07,  1.90e+07; ...
                +1,        +1];
    p.A_F_c2 = [ 5.01e-05; ...
                 2.05e+02; ...
                 3.03e+04; ...
                +1];
%     p.A_F_c1 = [ 5.4537e-1,  5.9384e-1; ...     % ha, hb
%                  1.2706e+2,  1.9524e+1; ...     % Ta, Tb
%                  2.8290e+0,  7.1057e-1];        % pa, pb
%     p.A_F_c2 = [ 1.2432e+0; ...                 % h
%                  1.8924e+2; ...                 % T
%                  1.1510e+0];                    % p
    
    % Other paramaters related to vesicle state changes
    p.rfr_T =   6.34;   % msec  - release refractory period
    
    p.dep_T = 1.00e+1;  % msec  - depression rate decay time constant
    p.h_dep = 0*4.94e-2;  % ms^-1 - depression rate increment
    p.k_d_0 = 0.00e+0;  % ms^-1 - equilibrium depression rate
    p.dsg_h = 2.00e+0;  % msec  - depression sigmoid half rise time
    p.dsg_p = 5.00e+0;  % depression sigmoid power
    
    p.act_T = 3.20e+2;  % msec  - reactivation rate decay time constant
    p.h_act = 7.55e-5;  % ms^-1 - reactivation rate increment
    p.k_a_0 = 5.00e-4;  % ms^-1 - equilibrium reactivation rate
    p.asg_h = 2.00e+0;  % msec  - activation sigmoid half rise time
    p.asg_p = 5.00e+0;  % activation sigmoid power
    
    p.rec_T = Inf;      % msec  - vesicle recycling time constant
    p.prm_T =   1.0e3;  % msec  - vesicle priming "fast" time constant
%     p.rec_T = 100e3;    % msec  - vesicle recycling time constant
%     p.prm_T =  23.6e3;  % msec  - vesicle priming "fast" time constant
    
    p.P_kar =   1;      % P(kiss-and-run)
    p.kar_T =   0.8e3;      % msec  - kiss-and-run time constant
%     p.P_kar =   0.44;   % P(kiss-and-run)
%     p.kar_T = 860;      % msec  - kiss-and-run time constant
    p.fil_T =  10.7e3;  % msec  - delayed vesicle retrieval time
    
    p.rdy_T =   1.0e3;      % msec  - vesicle priming "slow" time constant
    p.rtr_T = Inf;      % msec  - time constant for vesicle to axon
    p.ant_T = Inf;      % msec  - time constant for vesicle from axon
%     p.rdy_T = 105.6e3;  % msec  - vesicle priming "slow" time constant
%     p.rtr_T = 500e3;    % msec  - time constant for vesicle to axon
%     p.ant_T = 500e3;    % msec  - time constant for vesicle from axon
end

spks = spks(:);

% Initialize vesicle pools.
RRP = p.RRP0;   % RPP <-> Readily Releasable Pool -> KAR/EMT
KAR = p.KAR0;   % RRP -> Kiss-And-Run pool -> PMP
EMT = p.EMT0;   % RRP -> EMpTy pool -> PMP
PMP = p.PMP0;   % KAR -> PreMature Pool -> RPP
RPP = p.RPP0;   % PMP -> Readily Priming Pool <-> RRP

% Accomodate multiple active zones.
if isempty(nAZ)
    nAZ = 1;
end
RRP = RRP * ones(1, nAZ);

vp_init = struct('KAR', KAR, 'EMT', EMT, ...
    'PMP', PMP, 'RPP', RPP, 'RRP', RRP);

% Set up parameters for vesicle release profiles.

% ------------ Synchronous Parameters ------------ %

% Profile
SC_f = p.S_C;
ST_f = p.S_T;
% Sh_f = p.S_h;
Sp_f = p.S_p;

% Facilitation
SFc1 = [0 0 1]; SFc2 = [0 0 1]; SFc3 = [0 1];
SFt1 = [0 0 1]; SFp0 = [0 1];

% ------------ Asynchronous Parameters ------------ %

% Profile
AC_f = p.A_C;
% AT_f = p.A_T;
% Ah_f = p.A_h;
% Ap_f = p.A_p;

% Facilitation
AFc1 = [0 0 1]; AFc2 = [0 1];

% Depression
k_RID = p.k_d_0;    % release-independent depression rate
k_FDR = p.k_a_0;    % frequency-dependent recovery from depression rate

% ------------------------- Pool and Event Data ------------------------ %

% "#define" vesicle pool IDs (for event dependencies).
SPK_0 = 0;
RRP_V = 1; RRP_D = 2; RRP_R = 3; RRP_P = 4;
KAR_0 = 5; EMT_0 = 6; PMP_0 = 7; RPP_0 = 8;

% "#define" event IDs.
spk =  1;

rel =  2; rfr =  3;
dep =  4; act =  5;
rec =  6; prm =  7;

kar =  8; fil =  9; rdy =  10;
rtr = 11; ant = 12;

% Next event times for entire synapse
next = zeros(12, 1) + Inf;

% Spike times
lastSpk = -Inf;

% -------------------------- RRP-Specific Data ------------------------- %

% "#define" RRP state IDs.
V = 1;  % vacant
D = 2;  % depressed
R = 3;  % refracted
P = 4;  % primed

% Next event times for each RRP active zone
next_AZ = zeros(6, nAZ) + Inf;  % rel, rfr, dep, act, rec, prm

% The next active zone to be involved in an event
AZ_eIDs = zeros(6, 1);

% The next release mechanism at each active zone
rel_mech = zeros(1, nAZ);

% The next reuptake mechanism and release component at each active zone
reuptake = zeros(1, nAZ);
rel_comp = zeros(1, nAZ);

% ------------------- Keep track of all event times. ------------------- %

% rels = 0*spks; Prel = 0*spks;
history = zeros(5*(1+length(spks)), 6);     % time, event ID, active zone,
                                            % release mech, reuptake mech,
                                            % release component
ns = 0; nh = 0;

% ------------- Run event-driven simulation until the end. ------------- %

if isempty(T)
    if ~isempty(spks)
        T = spks(end) + spks(1);
    else
        T = 1000;
    end
end
t = 0; evt = 0; az = 1;     % current time, event, and active zone
while t < T
    if evt > 0
        nh = nh + 1;
        
        % Record for posterity.
        history(nh,1:2) = [t, evt];
        
        % Ensure that current event does not interfere.
        next(evt) = Inf;
        
        % Keep track of individual active zones.
        if any(evt == [rel, rfr, dep, act, rec, prm])
            az = AZ_eIDs(evt-1);    % active zone that has this event
            
            % Ensure that current event does not interfere.
            next_AZ(evt-1,az) = Inf;
            
            % Record for posterity.
            history(nh,3) = az;
            if evt == rel
                history(nh,4) = rel_mech(az);   % release mechanism
                history(nh,5) = reuptake(az);   % reuptake mechanism
                history(nh,6) = rel_comp(az);   % release component
            end
        end
    end
    
    % Act based on which event occurred: update the state
    % of the synapse, and prepare for future events.
    switch evt
        case spk    % spike time
            % Perform facilitations of vesicle release parameters
            % and of depletion-independent depression rate.
            facil_all(t - lastSpk);
            
            % Update spike sequence.
            lastSpk = t;
            
            % Update state of the synapse, and prepare for future events.
            delta_pool(SPK_0, [], 0);
            
        case rel    % release time
            % Choose reuptake mechanism of released vesicle.
            if reuptake(az)
                delta_pool(KAR_0, az, +1);
            else
                delta_pool(EMT_0, az, +1);
            end
            
            % Remove vesicle from RRP and send all
            % active vesicles into refractory period.
            delta_pool(RRP_V, az, +1);
            delta_pool(RRP_R, az,  RRP(P,az) - 1);
            delta_pool(RRP_P, az, -RRP(P,az));
            
        case rfr    % end refractory period
            % Activate a refracted vesicle.
            delta_pool(RRP_R, az, -1);
            delta_pool(RRP_P, az, +1);
            
        case dep    % depress a vesicle
            % Deactivate a vesicle from the releasable state.
            delta_pool(RRP_D, az, +1);
            delta_pool(RRP_P, az, -1);
            
        case act    % recover from depression (activate)
            % Return vesicle to default state.
            delta_pool(RRP_D, az, -1);
            delta_pool(RRP_P, az, +1);
            
        case rec    % RPP <- RRP (recycle)
            % Recycle a vesicle from the releasable state.
            delta_pool(RRP_V, az, +1);
            delta_pool(RRP_P, az, -1);
            delta_pool(RPP_0, az, +1);
            
        case prm    % RPP -> RRP (prime)
            % Prime a vesicle for release.
            delta_pool(RRP_V, az, -1);
            delta_pool(RRP_P, az, +1);
            delta_pool(RPP_0, az, -1);
            
        case kar    % KAR -> PMP (kiss-and-run)
            % Refill a kiss-and-run vesicle.
            delta_pool(KAR_0, [], -1);
            delta_pool(PMP_0, [], +1);
            
        case fil    % EMT -> PMP (fill)
            % Refill a fully fused vesicle.
            delta_pool(EMT_0, [], -1);
            delta_pool(PMP_0, [], +1);
            
        case rdy    % PMP -> RPP (ready)
            % Ready a filled vesicle for priming.
            delta_pool(PMP_0, [], -1);
            delta_pool(RPP_0, [], +1);
            
        case rtr    % AXN <- PMP (retrograde)
            % Remove a vesicle from the synapse to the axon.
            delta_pool(PMP_0, [], -1);
            
        case ant    % AXN -> PMP (anterograde)
            % Deliver a new vesicle to the synapse from the axon.
            delta_pool(PMP_0, [], +1);
            
        otherwise
            % Initialize next event times.
            delta_pool(SPK_0, 1:nAZ, 0);
            delta_pool(RRP_V, 1:nAZ, 0);
            delta_pool(RRP_D, 1:nAZ, 0);
            delta_pool(RRP_R, 1:nAZ, 0);
            delta_pool(RRP_P, 1:nAZ, 0);
            delta_pool(KAR_0, 1:nAZ, 0);
            delta_pool(EMT_0, 1:nAZ, 0);
            delta_pool(PMP_0, 1:nAZ, 0);
            delta_pool(RPP_0, 1:nAZ, 0);
    end
    
    % Find earliest RRP event times and which active zones they come from.
    [next(rel:prm), AZ_eIDs] = min(next_AZ, [], 2);
    
    % Choose earliest next event time, and jump to it.
    [t, evt] = min(next);
end

% Eliminate excess entries.
history = history(1:nh,:);

% Create final output.
rels = history(history(:,2)==rel,1);
Prel = 1 + 0*rels;

% -------------------------- Helper Functions -------------------------- %

% Define how pool sizes change with events and
% how that affects the timing of future events.
    function delta_pool(pool, AZid, delta)
        
        % Update pool size and gather dependencies.
        [events, pool_sz, AZs] = dependencies(pool, AZid, delta);
        
        % Update next event times for the new vesicles added to the pool.
        if delta > 0
            for e = events
                % Calculate next event time just for the new vesicle(s).
                % If they act before the old vesicles, use the new time.
                next_time(e, AZs, min(delta, pool_sz), t);
            end
        
        % Update next event times for the vesicles still left in the pool.
        elseif delta < 0
            for e = events
                if e == evt
                    % Event just completed, so calculate
                    % the next time for this event.
                    next_time(e, AZs, pool_sz, t);
                elseif rand(1) < -delta/(pool_sz-delta)
                    % Removed vesicle(s) (delta) would have done this event
                    % next, so calculate when remaining vesicles will act.
                    % ??? next(e) = Inf; because it was removed?
                    next_time(e, AZs, pool_sz, next(e));
                end
            end
            
        % No change to a vesicle pool size, just a change in rates.
        else
            for e = events
                next_time(e, AZs, pool_sz, t);
            end
        end
        
    end

% Define dependencies of events on vesicle pools.
    function [events, pool_sz, AZs] = dependencies(pool, AZid, delta)
        
        AZs = AZid;
        switch pool
            case SPK_0
                events = [spk, rel, dep];
                pool_sz = RRP(P,:);
                AZs = 1:nAZ;    % release and depress at all active zones
            case RRP_V
                events = prm;
                RRP(V,AZid) = RRP(V,AZid) + delta;
                pool_sz = min(RPP, RRP(V,AZid));
            case RRP_D
                events = act;
                RRP(D,AZid) = RRP(D,AZid) + delta;
                pool_sz = RRP(D,AZid);
            case RRP_R
                events = rfr;
                RRP(R,AZid) = RRP(R,AZid) + delta;
                pool_sz = RRP(R,AZid);
            case RRP_P
                events = [rel, dep, rec];
                RRP(P,AZid) = RRP(P,AZid) + delta;
                pool_sz = RRP(P,AZid);
            case KAR_0
                events = kar;
                KAR = KAR + delta;
                pool_sz = KAR;
            case EMT_0
                events = fil;
                EMT = EMT + delta;
                pool_sz = EMT;
            case PMP_0
                events = [rdy, rtr];
                PMP = PMP + delta;
                pool_sz = PMP;
            case RPP_0
                events = prm;
                RPP = RPP + delta;
                pool_sz = min(RPP, RRP(V,:));
                AZs = 1:nAZ;    % priming to all active zones
        end
        
    end

% Define how to choose subsequent event times for the whole synapse.
    function nt = next_time(event, AZs, pool_sz, time)
        
        is_RRP = false; nt = Inf;
        
        % Calculate next event time.
        if event == spk     % spike time
            if ns < length(spks)
                ns = ns + 1;
                nt = spks(ns);
            else
                nt = Inf;
            end
        elseif pool_sz > 0
            switch event
%                 case spk    % spike time
%                     if ns < length(spks)
%                         ns = ns + 1;
%                         nt = spks(ns);
%                     else
%                         nt = Inf;
%                     end
%                     ns
                    
                case {rel, rfr, dep, act, rec, prm}     % RRP events
                    is_RRP = true;
                    for z = 1:length(AZs)
                        if length(pool_sz) == length(AZs)
                            ps = pool_sz(z);
                        else
                            ps = pool_sz;
                        end
                        next_RRP_time(event, AZs(z), ps, time);
                    end
                    
                case kar    % KAR -> PMP (kiss-and-run)
                    nt = time + exprnd(p.kar_T / pool_sz);
                    
                case fil    % EMT -> PMP (fill)
                    nt = time + exprnd(p.fil_T / pool_sz);
                    
                case rdy    % PMP -> RPP (ready)
                    nt = time + exprnd(p.rdy_T / pool_sz);
                    
                case rtr    % AXN <- PMP (retrograde)
                    nt = time + exprnd(p.rtr_T / pool_sz);
                    
                case ant    % AXN -> PMP (anterograde)
                    nt = time + exprnd(p.ant_T);
            end
        end
        
        % Reassign next event time if it comes before current.
        if ~is_RRP
            if nt < next(event)
                next(event) = nt;
            end
        end
    end

% Define how to choose subsequent event times for the RRP.
    function next_RRP_time(event, AZid, pool_sz, time)
        
        nt = Inf;
        
        % Calculate next event time.
        if pool_sz > 0
            switch event
                case rel    % release time
                    [S_rel, S_cmp] = ...
                        Ca_evt_time(SC_f, ST_f, ...
                        p.S_h, Sp_f, p.sr0, pool_sz, time);
                    [A_rel, A_cmp] = ...
                        Ca_evt_time(AC_f, p.A_T, ...
                        p.A_h, p.A_p, p.ar0, pool_sz, time);
                    [nt, rm] = min([S_rel, A_rel]);
                    
                    % Track release only if it's first.
                    if nt < next_AZ(rel-1,AZid)
                        next_AZ(rel-1,AZid) = nt;
                        rel_mech(AZid) = rm;
                        reuptake(AZid) = (rand(1) < p.P_kar);
                        if rm == 1
                            rel_comp(AZid) = S_cmp;
                        else
                            rel_comp(AZid) = A_cmp;
                        end
                    end
                    
                case rfr    % end refractory period
                    nt = time + exprnd(p.rfr_T / pool_sz);
                    
                case dep    % depress a vesicle
                    nt = Ca_evt_time(k_RID, p.dep_T, ...
                        p.dsg_h, p.dsg_p, p.k_d_0, pool_sz, time);
                    
                case act    % recover from depression (activate)
                    nt = Ca_evt_time(k_FDR, p.act_T, ...
                        p.asg_h, p.asg_p, p.k_a_0, pool_sz, time);
                    
                case rec    % RPP <- RRP (recycle)
                    nt = time + exprnd(p.rec_T / pool_sz);
                    
                case prm    % RPP -> RRP (prime)
                    nt = time + exprnd(p.prm_T / pool_sz);
            end
        end
        
        % Reassign next event time if it comes before current.
        if event ~= rel && nt < next_AZ(event-1,AZid)
            next_AZ(event-1,AZid) = nt;
        end
        
    end

% Define vesicle release timing function.
    function [rt, cmp] = Ca_evt_time(C, T, h, p, r0, pool_sz, time)
        
        evt_C = time - lastSpk;
        while true
            % Find out when each component of this
            % release mechanism wants to release.
            evt_C = nXPrnd(pool_sz, C(:), T(:), evt_C);
            
            % Take the earliest of the resulting release times.
            [evt_C, cmp] = min(evt_C);
            
            % Determine whether the vesicle that would want to release at
            % these times actually sees the Ca from the current spike.
            P_new = 1 / (1 + (evt_C / h)^(-p));
            if rand(1) < P_new
                evt_C = lastSpk + evt_C;
                break
            end
        end
        
        % Calculate spike-independent release time.
        evt_M = time + exprnd(1 / (pool_sz * r0));
        
        % Return the earliest release time among choices.
        [rt, sm] = min([evt_C; evt_M]);
        if sm == 2
            cmp = 0;
        end
        
    end

% Define function to facilitate all profile parameters.
    function facil_all(del_T)
        
        % ------------ Synchronous Parameters ------------ %
        SFc1 = facil_param(del_T, SFc1, p.S_F_c1);
        SC_f(1) = p.S_C(1) * SFc1(end);
        
        SFt1 = facil_param(del_T, SFt1, p.S_F_t1);
        ST_f(1) = p.S_T(1) * SFt1(end);
        
        SFc2 = facil_param(del_T, SFc2, p.S_F_c2);
        SC_f(2) = p.S_C(2) * SFc2(end);
        
        SFc3 = facil_param(del_T, SFc3, p.S_F_c3);
        SC_f(3) = p.S_C(3) * SFc3(end);
        
        SFp0 = facil_param(del_T, SFp0, p.S_F_p0);
        Sp_f = p.S_p * SFp0(end);
        
        % ------------ Asynchronous Parameters ------------ %
        AFc1 = facil_param(del_T, AFc1, p.A_F_c1);
        AC_f(1) = p.A_C(1) * AFc1(end);
        
        AFc2 = facil_param(del_T, AFc2, p.A_F_c2);
        AC_f(2) = p.A_C(2) * AFc2(end);
        
        % Facilitate release-independent depression rate.
        k_RID = p.h_dep + k_RID * exp(-del_T / p.dep_T);
        
        % Facilitate frequency-dependent recovery from depression rate.
        k_FDR = p.h_act + k_FDR * exp(-del_T / p.act_T);
        
    end

% Define ISI-dependent facilitation function.
    function F_ratio = facil_param(del_T, F_ratio, F_vars)
        
        F_ratio(end) = 1;
        for j = 1:size(F_vars,2)
            phi_F = F_vars(1,j);
            tau_F = F_vars(2,j);
            max_F = F_vars(3,j);
            pow_F = F_vars(4,j);
            
            % f(n+1) = (f(n) + phi*(1 - f(n))) * exp(-ISI / T)
            F_ratio(j) = (F_ratio(j) + phi_F * (1 - F_ratio(j)));
            F_ratio(j) = F_ratio(j) .* exp(-del_T ./ tau_F);
            
            % F(n+1) = (f(n+1) + 1)^X
            F_ratio(end) = F_ratio(end) * (max_F * F_ratio(j) + 1)^pow_F;
        end
        
    end

end

%  % Ca-dependent priming rate profile parameters
%  p.pr0 = 1.0000e-1;  %  uM   - steady-state Ca concentration
%  p.P_C = ...         % Ca-spike profile coefficients
%      [8.6947e+1; ...     % fast (uM)
%       8.0572e+0];        % slow (uM*msec)
%  p.P_T = ...         % msec  - double exponential time constants
%      [6.8204e-1; ...     % fast (msec)
%       2.8788e-2];        % slow (msec^-1)
%  p.P_h = 1.0419e+0;  % msec  - initial ramp half-rise time
%  p.P_p = 7.7942e+0;  % power on sigmoid
%  
%  % Ca-dependent priming rate facilitation parameters
%  p.P_F_c1 = [ 1.2974e-2; ...                 % h
%               3.1327e+1; ...                 % T
%               1.3339e+1];                    % p
%  p.P_F_c2 = [ 6.5512e-1; ...                 % h
%               1.4716e+2; ...                 % T
%               1.4127e+0];                    % p
%  
%  % release parameters
%  p.T_refr
%  
%  p.k_pr = 0.5e-3;  % msec^-1   - vesicle prep rate
% 
%  p.k_d0 = 0.1e-3;  % msec^-1   - base vesicle priming rate
%  p.k_dC = 1.0e-3;  % 1/(uM*ms) - Ca-dependent vesicle priming rate
%  p.k_ud = (N_RRP_m - RRP0) * (k_d0 + k_dC * Ca_0) / RRP0;  % undock rate
%  
%  % vesicle retrieval: kiss-and-run, compensatory, evoked
%  p.T_karn
%  p.T_comp
%  p.comp_M = log(T_comp / sqrt(72 * 350^2 / T_comp^2 + 1));
%  p.comp_S = sqrt(log(72 * 350^2 / T_comp^2 + 1));
%  p.T_evok =   424; % msec - spike-evoked vesicle retrieval
%  p.fr_kce = [0.44; 0.20; 0.36];
% 
% % Find out when the previous spike would cause release, if any.
% if 1 %evt ~= spk  % Wait for release to use this.
%     % Find out when each component of this
%     % release mechanism wants to release.
%     rel_P = nXPrnd(n, Cp(:), Tp(:), time - prevSpk);
%     
%     % Take the earliest of the resulting release times.
%     [rel_P, cmp_P] = min(rel_P);
%     prev_R = prevSpk + rel_P;
% else
%     prev_R = Inf;
%     cmp_P = 0;
% end