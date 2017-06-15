% Jonathan Garcia
% 2016-12-26

classdef (ConstructOnLoad) EventTime < event.EventData
    properties
        new_time    % time for new event to be added to Scheduler
    end
    
    methods
        % Constructor
        function ET = EventTime(time)
            ET.new_time = time;
        end
    end
end