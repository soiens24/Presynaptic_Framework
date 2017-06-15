% Jonathan Garcia
% 2016-12-25

classdef VesiclePool < handle
    properties
        pool_size   % number of vesicles currently in the pool
        capacity    % maximum number of vesicles allowed in pool
        shared      % number of vesicles currently in another state
    end
    
    events
        PoolUpdate  % vesicle pool has changed in size
    end
    
    methods
        % Constructor
        function vph = VesiclePool(sz, mx, sh)
            if nargin < 3
                sh = 0;     % any vesicles in current state
            if nargin < 2
                mx = Inf;   % no specified upper bound
            if nargin < 1
                sz = 0;     % empty by default
            end
            end
            end
            vph.pool_size = min(sz, mx);
            vph.capacity = max(mx, sz);
            vph.shared = min(sh, mx-sz);
        end
        
        % When the pool changes size, notify all processes that
        % listen to it, which use this pool as a vescile source.
        function update(vph, delta, type, new_time, share_TF)
            if delta <= 0
                original = vph.pool_size;   % vesicles
            else
                original = vph.vacancies;   % vacancies
            end
            vph.pool_size = vph.pool_size + delta;
            if nargin < 5
                share_TF = false;
            end
            if share_TF
                % Do this if this update is shared with another pool
                % (i.e., the vesicles are not really moving, just changing
                % state within the same number of limited sites).
                vph.shared = vph.shared - delta;
            end
            notify(vph, 'PoolUpdate', ...
                DeltaPool(delta, type, new_time, original, share_TF));
        end
        
        % Return the number of vacancies open for new vesicles.
        function N_vac = vacancies(vph, share_TF)
            if nargin < 2
                share_TF = false;
            end
            if share_TF
                N_vac = vph.capacity - vph.pool_size;
            else
                N_vac = vph.capacity - (vph.pool_size + vph.shared);
            end
        end
    end
end