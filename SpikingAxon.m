% Jonathan Garcia
% 2016-12-26

classdef SpikingAxon < handle
    properties
        % AZ_list     % list of all ActiveZone's that come off this axon
        % AZ_delays   % list of delays from spike onset to spike arrival
        lastSpk     % timing of the most recent spike
        del_T       % interspike interval that preceded most recent spike
    end
    
    events
        Spike
    end
    
    methods
        % Constructor
        function sah = SpikingAxon()
            sah.lastSpk = -Inf; % There was no last spike.
            sah.del_T = Inf;    % Synapse has been at rest forever.
        end
        
        % Generate a spike at a scheduled time.
        function spike(sah, new_time)
            sah.del_T = new_time - sah.lastSpk;
            sah.lastSpk = new_time;
            notify(sah, 'Spike', EventTime(new_time));
        end
    end
end