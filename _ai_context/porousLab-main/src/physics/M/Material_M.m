%% Material_M Class
% This class represents a material model that combines porous media 
% properties with mechanical constitutive behavior. It provides methods to 
% evaluate the mechanical constitutive law, retrieve the number of state 
% variables, and check for plastic strain behavior.
%
%% Methods
% * *mechanicalLaw*: Evaluates the mechanical constitutive law at a given 
%                    integration point to return the stress tensor and the 
%                    material stiffness matrix.
% * *getNumberStateVar*: Returns the number of state variables associated 
%                        with the mechanical constitutive law.
% * *hasPlasticStrain*: Checks if the material exhibits elasto-plastic 
%                       behavior.
% 
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00: Initial version (January 2024).
%
%% Class definition
classdef Material_M < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        porousMedia = PorousMedia();
        mechanical  = [];
    end  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Material_M(matData,lc)
            this.porousMedia = matData.porousMedia;
            % Mechanical constitutive behavior
            if strcmp('elastic',matData.porousMedia.mechanical)
                this.mechanical = MechanicalLinearElastic();
            elseif strcmp('vonMises',matData.porousMedia.mechanical)
                this.mechanical = MechanicalElastoPlasticVonMises();
            elseif strcmp('druckerPrager',matData.porousMedia.mechanical)
                this.mechanical = MechanicalElastoPlasticDruckerPrager();
            elseif strcmp('elasticDruckerPrager',matData.porousMedia.mechanical)
                this.mechanical = MechanicalNonlinearElasticDruckerPrager();
            elseif strcmp('mohrCoulomb', matData.porousMedia.mechanical)
                this.mechanical = MechanicalElastoPlasticMohrCoulomb();
            elseif strcmp('nonlinearAsymptotic',matData.porousMedia.mechanical)
                this.mechanical = MechanicalNonlinearAsymptotic();
            elseif strcmp('isoDamage',matData.porousMedia.mechanical)
                this.mechanical = MechanicalIsotropicDamage(lc);
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
        % Checks if the material exhibits elasto-plastic behaviour
        function flag = hasPlasticStrain(this)
            flag = this.mechanical.isElastoPlastic();
        end

    end
end