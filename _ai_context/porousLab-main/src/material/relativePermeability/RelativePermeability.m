%% RelativePermeability Class
% This class defines an abstract base class for modeling relative 
% permeability in porous media. It provides a framework for implementing 
% specific relative permeability models by defining an abstract method 
% _calculate_ that must be implemented in derived classes.
%
%% Methods
% *calculate*: Abstract method to compute the relative permeability. This 
%              method must be implemented in subclasses.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class Definition
classdef RelativePermeability < handle    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        id = 'name1';   
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RelativePermeability(matModel)
            this.id = matModel;
        end
    end

    %% Abstract methods
    methods(Abstract)

        %------------------------------------------------------------------
        % Compute the relative permeability
        kr = calculate(this, Sl, porousMedia);

        %------------------------------------------------------------------
        % Compute the derivative of the relative permeability wrt the
        % saturation
        dkrdSl = derivative(this, Sl, porousMedia);
        
    end

end