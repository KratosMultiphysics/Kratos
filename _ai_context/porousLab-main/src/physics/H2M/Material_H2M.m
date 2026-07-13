%% Material_H2M Class
% This class extends the `Material_H2` class to include mechanical 
% constitutive behavior for porous media. It provides methods to evaluate 
% mechanical laws, retrieve state variables, and compute mechanical 
% compressibility coefficients.
%
%% Methods
% * *mechanicalLaw*: Evaluates the mechanical constitutive law at a given 
%                    integration point and returns the stress and the 
%                    constitutive matrix.
% * *getNumberStateVar*: Returns the number of state variables associated 
%                        with the mechanical constitutive law.
% * *biotCoeff*: Returns the Biot coefficient of the porous media.
% * *mechanicalCompressibilityCoeffs*: Computes the mechanical 
%                                      compressibility coefficients for 
%                                      liquid and gas phases based on the 
%                                      saturation `Sl`.
% * *hasPlasticStrain*: Checks if the material exhibits elasto-plastic 
%                       behavior.
% 
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class definition
classdef Material_H2M < Material_H2    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        mechanical = [];
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Material_H2M(matData)
            this = this@Material_H2(matData);
            % Mechanical constitutive behavior
            if strcmp('elastic',matData.porousMedia.mechanical)
                this.mechanical = MechanicalLinearElastic();
            end
        end
    end
    %% Public methods
    methods
        % -----------------------------------------------------------------
        % Evaluate the mechanical constitutive law
        function [stress,D] = mechanicalLaw(this,ip)
            [stress,D] = this.mechanical.eval(this.porousMedia,ip);
        end

        % -----------------------------------------------------------------
        % Get the number of state variables associated with the mechanical
        % constitutive law
        function nstVar = getNumberStateVar(this)
            nstVar = this.mechanical.nstVar;
        end

        % -----------------------------------------------------------------
        % Returns the biot coefficient
        function biot = biotCoeff(this)
            biot = this.porousMedia.biot;
        end

        % -----------------------------------------------------------------
        % Evaluate the mechanical constitutive law
        function [cul,cug] = mechanicalCompressibilityCoeffs(this,Sl)
            cul = this.porousMedia.biot * Sl;
            cug = this.porousMedia.biot * (1.0 - Sl);
        end

        % -----------------------------------------------------------------
        % Checks if the material exhibits elasto-plastic behaviour
        function flag = hasPlasticStrain(this)
            flag = this.mechanical.isElastoPlastic();
        end
    end
end