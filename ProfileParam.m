% Jonathan Garcia
% 2016-12-25

classdef ProfileParam
    properties
        % baseline parameter value
        param_0     % any of: r0, hh, pp, C1, T1, ... Cn, Tn
        
        % facilitation parameters
        fs      % internal facilitating parameter
        Fn      % facilitation factor
        facil_TNX   % each row: decay time constant (T),
                    %           steps to saturation (N),
                    %           non-linear exponent (X).
    end
    
    methods
        % Constructor
        function pph = ProfileParam(p0, Ts, Ns, Xs)
            pph.param_0 = p0;
            if nargin > 1   % parameter facilitates
                if numel(T) ~= numel(N) || numel(N) ~= numel(X)
                    error('T, N, and X must have same size.');
                end
                pph.facil_TNX = [Ts(:), Ns(:), Xs(:)];
                pph.fs = 0*Ts(:);
                pph.Fn = 1;
            else
                pph.facil_TNX = []; % no facilitation
                pph.fs = [];
            end
            pph.Fn = 1;     % Baseline is default.
        end
        
        % Perform facilitation.
        function pph = facilitate(pph, del_T)
            if isempty(pph.facil_TNX)
                return  % nothing to do
            end
            
            % Go through each facilitation component individually.
            pph.Fn = 1;
            for j = 1:length(pph.fs)
                T = pph.facil_TNX(j,1);     % decay time constant (T)
                N = pph.facil_TNX(j,2);     % steps to saturation (N)
                X = pph.facil_TNX(j,3);     % non-linear exponent (X)
                
                g = pph.fs(j) * exp(-del_T / T);    % left-overs
                pph.fs(j) = g + 1 - (g/N)^N;        % limited facilitation
                
                pph.Fn = pph.Fn * pph.fs(j)^X;      % combining components
            end
        end
        
        % Access the current value of the parameter.
        function pv = value(pph)
            % product of baseline value and facilitation factor
            pv = pph.param_0 * pph.Fn;
        end
    end
end