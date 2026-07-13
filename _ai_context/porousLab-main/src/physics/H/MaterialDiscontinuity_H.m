%% MaterialDiscontinuity_H class
% This class represents a material discontinuity in a porous medium, 
% characterized by its initial aperture and the fluid properties. 
% It provides methods to compute the longitudinal permeability 
% coefficient and compressibility based on the material's properties.
%
%% Methods
% * *longitudinalPermeability*: Computes the longitudinal permeability 
%                               coefficient based on the cubic law.
% * *compressibility*: Computes  the compressibility of the material 
%                      discontinuity based on its aperture and fluid 
%                      properties.
% 
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class definition
classdef MaterialDiscontinuity_H < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        initialAperture = 0.0;
        leakoff         = 1.0;
        liquidFluid     = Fluid();
    end  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MaterialDiscontinuity_H(matData)
            this.initialAperture = matData.initialAperture;
            this.liquidFluid = matData.liquidFluid;
            this.leakoff = matData.leakoff;
        end
    end
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Computes the longitudinal permeability coefficient
        function kl = longitudinalPermeability(this)
            w = this.initialAperture();
            kl = this.cubicLaw(w,this.liquidFluid.mu);
        end

        %------------------------------------------------------------------
        % Computes  the compressibility of the material discontinuity 
        % based on its aperture and fluid properties.
        function c = compressibility(this)
            w = this.initialAperture();
            c = w / this.liquidFluid.K;
        end

    end
    %% Static methods
    methods(Static)
        %------------------------------------------------------------------
        % Computes the longitudinal permeability coefficient based on the
        % cubic's law
        function kl = cubicLaw(w,mu)
            kl = w*w*w/12.0/mu;
        end
    end
end