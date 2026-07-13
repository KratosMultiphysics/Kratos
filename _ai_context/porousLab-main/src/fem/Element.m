%% Element Class
% This is an abstract class that defines a finite element.
%
%% Authors
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
% 
%% Class definition
classdef Element < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        type = [];
    end

    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Element(type)
            if nargin > 0
                this.type = type;
            end
        end
    end
end
