classdef Dendrite
    properties
        ID                % Unique identifier
        params            % Parameters: alpha, b, Tau, TauR, NaX, gc
        InitialConditions % Initial conditions [u0; v0]
    end
    
    methods
        function obj = Dendrite(ID, params, InitialConditions)
            if nargin < 1 || isempty(ID)
                error('Dendrite ID must be specified.');
            else
                obj.ID = ID;
            end
            if nargin < 2 || isempty(params)
                % Default parameters
                obj.params.alpha = 10;
                obj.params.b = 2;
                obj.params.Tau = 1700;
                obj.params.TauR = 1700;
                obj.params.NaX = 0.05;
                obj.params.gc = 4;
            else
                obj.params = params;
            end
            if nargin < 3 || isempty(InitialConditions)
                % Default initial conditions
                obj.InitialConditions = [0.2; 0.2];
            else
                obj.InitialConditions = InitialConditions;
            end
        end
    end
end
