% Jonathan Garcia
% 2016-12-25

classdef (ConstructOnLoad) DeltaPool < event.EventData
    properties
        delta   % amount by which the vesicle pool changes size
        type    % SynProcess handle for the event type which caused delta
        new_time    % time at which vesicle pool changes size
        P_taken % probability that a removed vesicle would have
                % participated in another process = # removed / # original
        share_TF    % whether this update is only a change of state, 
                    % rather than a change of location, affecting vacancy
    end
    
    methods
        % Constructor
        function Dpool = DeltaPool(delta, type, new_time, original, share_TF)
            Dpool.delta = delta;
            Dpool.type  = type;
            Dpool.new_time = new_time;
            Dpool.P_taken = abs(delta) / original;
            Dpool.share_TF = share_TF;
        end
    end
end