%% MaterialDiscontinuity_M Class
% This class represents a material discontinuity with mechanical and 
% constitutive behavior. It provides methods to evaluate the mechanical 
% constitutive law, retrieve the number of state variables, and check 
% for plastic strain behavior.
%
%% Methods
% * *mechanicalLaw*: Evaluates the mechanical constitutive law at a given 
%                    integration point and returns the stress and material 
%                    stiffness matrix.
% * *getNumberStateVar*: Returns the number of state variables associated 
%                        with the mechanical constitutive law.
% * *hasPlasticStrain*: Checks if the material exhibits elasto-plastic 
%                       behavior.
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00: Initial version (January 2024).
%
%% Class definition
classdef MaterialDiscontinuity_M < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        parameters  = [];
        mechanical  = [];
    end  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MaterialDiscontinuity_M(matData)
            % Create the material struct
            this.parameters = struct( ...
                'initialAperture',    matData.initialAperture, ...
                'normalStiffness',    matData.normalStiffness, ...
                'shearStiffness',     matData.shearStiffness,...
                'cohesion',           matData.cohesion,...
                'frictionAngle',      matData.frictionAngle,...
                'dilationAngle',      matData.dilationAngle,...
                'tensionCutOff',      matData.tensionCutOff,...
                'contactPenalization',matData.contactPenalization,...
                'maximumClosure',     matData.maximumClosure);
            % Mechanical constitutive behavior
            if strcmp('elastic',matData.cohesiveLaw)
                this.mechanical = MechanicalCohesiveLinearElastic();
            elseif strcmp('mohrCoulomb',matData.cohesiveLaw)
                this.mechanical = MechanicalCohesiveMohrCoulomb();
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

    end
end