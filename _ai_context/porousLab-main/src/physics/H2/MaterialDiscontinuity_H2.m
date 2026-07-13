%% MaterialDiscontinuity_H2 Class
% This class represents a two-phase flow material discontinuity in a porous
% medium, characterized by its initial aperture, porosity, leak-off, and
% liquid/gas fluid properties. It extends _Material_H2_ and provides the
% longitudinal permeability of the discontinuity.
%
%% Methods
% * *longitudinalPermeability*: Computes the longitudinal permeability
%                               coefficient based on the cubic law.
% 
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class definition
classdef MaterialDiscontinuity_H2 < Material_H2
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        initialAperture = 0.0;
        leakoff         = 1.0;
        porosity        = 1.0;
    end  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MaterialDiscontinuity_H2(matData)
            this = this@Material_H2(matData)
            if isempty(matData.initialAperture) == false
                this.initialAperture = matData.initialAperture;
            end
            if isempty(matData.leakoff) == false
                this.leakoff = matData.leakoff;
            end
            if isempty(matData.porosity) == false
                this.porosity = matData.porosity;
            end
        end
    end
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Computes the longitudinal permeability coefficients based on the
        % cubic's law
        function k = longitudinalPermeability(this)
            w = this.initialAperture();
            k = this.cubicLaw(w);
        end

    end
    %% Static methods
    methods(Static)
        %------------------------------------------------------------------
        % Computes the longitudinal permeability coefficient based on the
        % cubic's law
        function kl = cubicLaw(w)
            kl = w*w*w/12.0;
        end
    end
end
