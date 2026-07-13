%% CapillaryPressure Class
% This class defines an abstract base class for modeling capillary pressure 
% and its associated properties in porous media. It provides a framework 
% for implementing specific capillary pressure-saturation relationships 
% and their derivatives.
%
%% Methods
% * *saturationDegree*: Computes the liquid saturation degree based on the 
%                       capillary pressure and porous media properties.
% * *derivativeSaturationDegree*: Computes the derivative of the liquid 
%                                 saturation degree with respect to the 
%                                 capillary pressure.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class Definition
classdef CapillaryPressure < handle    
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        id = 'name1';   
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = CapillaryPressure(matModel)
            this.id = matModel;
        end
    end

    %% Abstract methods
    methods(Abstract)

        %------------------------------------------------------------------
        % Liquid saturation degree
        Sw = saturationDegree(this, pc, porousMedia);

        %------------------------------------------------------------------
        % Derivative of the liquid saturation degree wrt to the capillary
        % pressure
        dSwdPc = derivativeSaturationDegree(this, pc, porousMedia);
        
    end

end