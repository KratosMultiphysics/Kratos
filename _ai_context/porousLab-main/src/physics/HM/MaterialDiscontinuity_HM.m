%% MaterialDiscontinuity_HM Class
% This class represents a hydro-mechanical material discontinuity,
% characterized by cohesive mechanical properties, aperture, leak-off, and
% liquid-fluid properties. It provides methods to evaluate the cohesive law,
% update aperture, and compute longitudinal permeability and
% compressibility.
%
%% Methods
% * *mechanicalLaw*: Evaluates the cohesive mechanical law.
% * *initializeAperture*: Initializes the aperture state variable.
% * *updateAperture*: Updates the aperture state variable.
% * *longitudinalPermeability*: Computes the longitudinal permeability
%                               coefficient based on the cubic law.
% * *compressibility*: Computes the compressibility of the material
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
classdef MaterialDiscontinuity_HM < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        initialAperture = 0.0;
        leakoff         = 1.0;
        liquidFluid     = Fluid();
        parameters      = [];
        mechanical      = [];
    end  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MaterialDiscontinuity_HM(matData)
            % Create the material struct
            this.parameters = struct( ...
                'initialAperture',    matData.initialAperture, ...
                'normalStiffness',    matData.normalStiffness, ...
                'shearStiffness',     matData.shearStiffness,...
                'contactPenalization',matData.contactPenalization);
            this.initialAperture = matData.initialAperture;
            this.liquidFluid = matData.liquidFluid;
            this.leakoff = matData.leakoff;
            % Mechanical constitutive behavior
            if strcmp('elastic',matData.cohesiveLaw)
                this.mechanical = MechanicalCohesiveLinearElastic();
            end
        end
    end
    %% Public methods
    methods

        % -----------------------------------------------------------------
        % Evaluate the mechanical constitutive law
        function [stress,D] = mechanicalLaw(this,ip)
            [stress,D] = this.mechanical.eval(this.parameters,ip);
        end

        % -----------------------------------------------------------------
        % Initialize the aperture state variable
        function initializeAperture(this,ip)
           this.mechanical.initializeAperture(this.parameters,ip);
        end

        % -----------------------------------------------------------------
        % Get the number of state variables associated with the mechanical
        % constitutive law
        function nstVar = getNumberStateVar(this)
            nstVar = this.mechanical.nstVar;
        end

        % -----------------------------------------------------------------
        % Check if the material is elasto-plastic or not
        function flag = hasPlasticStrain(this)
            flag = this.mechanical.isElastoPlastic();
        end

        % -----------------------------------------------------------------
        % Get the aperture from the previous state variables
        function w = getAperture(~,ip)
            w = ip.statevarOld(1);
        end

        % -----------------------------------------------------------------
        % Update the aperture state variable from the normal strain increment
        function updateAperture(~,ip)
            ip.statevar(1) = ip.statevarOld(1) + (ip.strain(2) - ip.strainOld(2));
        end

        %------------------------------------------------------------------
        % Computes the longitudinal permeability coefficient based on the
        % cubic's law
        function kl = longitudinalPermeability(this,ip)
            w = this.getAperture(ip);
            kl = w*w*w/12.0/this.liquidFluid.mu;
        end

        %------------------------------------------------------------------
        % Computes  the compressibility of the material discontinuity 
        % based on its aperture and fluid properties.
        function c = compressibility(this)
            w = this.initialAperture();
            c = w / this.liquidFluid.K;
        end

    end
end
