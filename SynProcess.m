% Jonathan Garcia
% 2016-12-25

classdef SynProcess < handle
    properties
        lastEvt     % most recent previous event time
        lastCmp     % most recent previous event component
        nextEvt     % next planned event time
        nextCmp     % next planned event component
        
        r0      % ms^-1 - spontaneous rate
        hh      % msec  - transition half-rise time
        pp      %   #   - transition steepness order
        Cs      % ms^-1 - exponential rate coefficients
        Ts      % msec  - exponential rate decay time constants
    end
    
    events
        EventUpdate     % This process has a new event time
                        % for the scheduler to handle.
        DecrementSource     % Take vesicle from source pool.
        IncrementTarget     % Donate vesicle to target pool.
    end
    
    methods
        % Constructor
        function sph = SynProcess(r0, hh, pp, Cs, Ts)
            sph.nextEvt = Inf;
            sph.set_params('r0', r0);
            if nargin > 1
                sph.set_params('hh', hh);
                sph.set_params('pp', pp);
                sph.set_params('Cs', Cs);
                sph.set_params('Ts', Ts);
            else
                sph.hh = ProfileParam(0);
                sph.pp = ProfileParam(0);
                sph.Cs = ProfileParam(0);
                sph.Ts = ProfileParam(0);
            end
        end
        
        % Make sure to construct parameters as ProfileParam's, not numbers.
        function set_params(sph, param_name, param_val)
            if isa(param_val, 'ProfileParam')
                sph.(param_name) = param_val;
            else
                sph.(param_name) = ProfileParam(param_val(1));
                for p = 2:length(param_val)
                    sph.(param_name)(p) = ProfileParam(param_val(p));
                end
            end
        end
        
        % Manually reset a parameter's facilitation function.
        function set_facils(sph, param_name, idx, T, N, X)
            if isempty(idx)
                idx = 1;
            end
            if numel(T) ~= numel(N) || numel(N) ~= numel(X)
                error('T, N, and X must have same size.');
            end
            sph.(param_name)(idx).facil_TNX = [T(:), N(:), X(:)];
            sph.(param_name)(idx).fs = 0*T(:);
        end
        
        % Get this process to respond to spike events from some source.
        function respond2spike(sph, spike_src, vesicle_src, event)
            % First, facilitate parameters relative to the last spike.
            sph.facilitate(spike_src.del_T);
            
            % Now, update next event after the current spike.
            sph.new_evt_time(vesicle_src.pool_size, ...
                event.new_time, spike_src.lastSpk);
        end
        
        % Get this process to handle updates to vesicle source pool size.
        function respond2pool(sph, spike_src, vesicle_src, event)
            % If the event passed was for this process, update
            % next event time regardless of the pool change.
            if event.type == sph
                sph.lastEvt = sph.nextEvt;
                sph.lastCmp = sph.nextCmp;
                sph.nextEvt = Inf;  % Ensure new event time is used.
                sph.new_evt_time(vesicle_src.pool_size, ...
                    event.new_time, spike_src.lastSpk);
                return
            end
            
            % If you're here, a process other than this one
            % changed the size of the vesicle source pool.
            
            % Update next event time for new vesicles added to source pool.
            if event.delta > 0
                % Calculate next event time just for the new vesicle(s).
                % If they act before the old vesicles, use the new time.
                sph.new_evt_time(event.delta, ...
                    event.new_time, spike_src.lastSpk);
                
            % Update next event time for vesicles left in source pool.
            elseif event.delta < 0
                % Update if the vesicle that would have participated next
                % in this process was taken away by another process.
                if rand(1) < event.P_taken
                    % Removed vesicle(s) (delta) would have done this event
                    % next, so calculate when remaining vesicles will act.
                    sph.nextEvt = Inf;    % Ensure new event time is used.
                    sph.new_evt_time(vesicle_src.pool_size, ...
                        event.new_time, spike_src.lastSpk);
                end
            end
        end
        
        % Perform facilitation for all parameters.
        function facilitate(sph, del_T)
            % Nothing happens here for most processes.
            sph.r0 = sph.r0.facilitate(del_T);
            sph.hh = sph.hh.facilitate(del_T);
            sph.pp = sph.pp.facilitate(del_T);  % facilitated (down) on p_S
            
            % Facilitate all exponential components.
            for j = 1:length(sph.Cs)
                sph.Cs(j) = sph.Cs(j).facilitate(del_T);
                sph.Ts(j) = sph.Ts(j).facilitate(del_T);
            end
            
            % % debugging
            % sph.current_state(true);
        end
        
        % Just checking what the rates are.
        function [Cs, Ts] = current_state(sph, disp_TF)
            Cs = zeros(size(sph.Cs));
            Ts = zeros(size(sph.Ts));
            for j = 1:length(Cs)
                Cs(j) = sph.Cs(j).value();
                Ts(j) = sph.Ts(j).value();
            end
            if disp_TF
                for j = 1:length(Cs)
                    fprintf('%.3e kHz\t', Cs(j));
                end
                fprintf('\n');
                for j = 1:length(Ts)
                    fprintf('%.3e ms\t', Ts(j));
                end
                fprintf('\n');
            end
        end
        
        % Update the next event time for this process if the new time comes
        % sooner than the previously calculated next event time.
        function new_evt_time(sph, pool_sz, time, lastSpk)
            [nt, cmp] = sph.evt_time(pool_sz, time, lastSpk);
            if nt < sph.nextEvt || sph.nextEvt == Inf
                sph.nextEvt = nt;
                sph.nextCmp = cmp;
                notify(sph, 'EventUpdate', EventTime(nt));
            end
        end
        
        % Calculate the next event time relative to pool size,
        % to the current time, and to the last spike time.
        function [nt, cmp] = evt_time(sph, pool_sz, time, lastSpk)
            nt = Inf;
            cmp = 0;
            if pool_sz == 0
                % No events possible with no vesicles in source pool.
                return
            end
            
            del_T = time - lastSpk;
            p = sph.pp.value(); 
            h_p = sph.hh.value()^p;
            while true
                % Find out when each component of this process wants to
                % act, and take the earliest of the resulting event times.
                evt_C = Inf;
                for j = 1:length(sph.Cs)
                    evt_J = SynProcess.nXPrnd(pool_sz, ...
                        sph.Cs(j).value(), ...
                        sph.Ts(j).value(), del_T);
                    if evt_J < evt_C
                        evt_C = evt_J;
                        cmp = j;
                    end
                end
                if evt_C == Inf
                    break   % no spike-evoked event times found
                end
                
                % Determine whether the vesicle that would want to act at
                % these times actually sees Ca from the current spike.
                t_p = evt_C^p;
                P_new = t_p / (t_p + h_p);
                if rand(1) < P_new
                    % The new event time is final. Move on.
                    evt_C = lastSpk + evt_C;
                    break
                else
                    % The new event time must come after this point, then.
                    del_T = evt_C;
                end
            end
            
            % Calculate spike-independent event time.
            evt_M = time + exprnd(1 / (pool_sz * sph.r0.value()));
            
            % Return the earliest event time among choices.
            [nt, sm] = min([evt_M; evt_C]);
            if sm == 1
                cmp = 0;
            end
        end
    end
    methods (Static)
        % Compute the next event time for a time-varying Poisson process
        % whose rate parameter decays exponentially with time.
        function next = nXPrnd(N, C, T, last)
            if nargin < 4
                last = 0;
            end
            
            % Nothing can happen before the exponential is set off.
            if last < 0
                last = 0;
            end
            
            % Determine whether another event occurs, and if so, when.
            val = exprnd(1 / (N * C * T * exp(-last / T)));
            
            if val >= 1
                next = Inf;     % no more events
            else
                next = last - T * log(1 - val);
            end
        end
    end
end