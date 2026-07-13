%% Anl Class
% This is an abstract class that defines an analysis object.
%
%% Authors
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
%
%% Class definition
classdef Anl < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        type = [];
    end

    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Anl(type)
            if nargin > 0
                this.type = type;
            end
        end
    end
end
