%% MechanicalNonlinearAsymptotic Class
% This class implements a nonlinear asymptotic mechanical law for material
% behavior. It extends the _MechanicalLinearElastic_ class and provides
% methods to compute the stress vector, constitutive matrix, shear modulus,
% and its gradient based on the material properties and strain invariants.
% 
% Reference:
% Pasquali, Paulo Roberto Zanella.
% "Análise limite de estruturas através de uma formulação em elasticidade
% não-linear." (2008).
%
%% Methods
% * *eval*: Computes the stress vector and the constitutive matrix for the 
%           material based on the material properties and strain
%           invariants.
% * *getShearModulus*: Computes the shear modulus based on the material 
%                      properties, volumetric strain, norm of the 
%                      deviatoric strain, and bulk modulus.
% * *getGradientShearModulus*: Computes the gradient of the shear modulus 
%                              with respect to the strain.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class Definition
classdef MechanicalNonlinearAsymptotic < MechanicalLinearElastic  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MechanicalNonlinearAsymptotic()
            this = this@MechanicalLinearElastic();
        end
    end

    %% Public methods
    methods
        %------------------------------------------------------------------
        % Compute the stress vector and the constitutive matrix
        function [stress,De] = eval(this,material,ip)

            % Decompose the strain tensor
            ev = this.volumetricStrain(ip.strain);
            ed = this.deviatoricStrain(ip.strain);

            % Norm deviatoric strain tensor
            ned = this.normDeviatoricStrain(ip.strain);

            % Material constants
            K = this.bulkModulus(material);

            % Shear modulus
            mu = this.getShearModulus(material,ev,ned,K);

            % Identity tensor
            iv = this.gradientStrainInvariantI1();
            Id = this.fourthOrderSymTensor();

            % Stress vector
            stress = K * ev * iv + 2.0 * mu * ed;

            % Shear modulus
            dmude = this.getGradientShearModulus(material,ip.strain,ev,ned,K);

            % Constitutive matrix
            De = K * (iv * iv');
            De = De + 2.0 * mu * (Id - 1.0/3.0 * (iv * iv'));
            De = De + 2.0 * dmude * (Id * ed)';
        end

        %------------------------------------------------------------------
        % Compute the shear modulus
        function mu = getShearModulus(~,material,ev,ned,K)
            if strcmp(material.asympt,'vonMises')
                mu = sqrt(2.0)*material.tauy/2.0;
            elseif  strcmp(material.asympt,'druckerPrager')
                mu = material.friction * (material.sy - K * ev) / 2.0;
            end
            mu = mu / (material.eref + ned);
        end

        %------------------------------------------------------------------
        % Compute the gradient of the shear modulus
        function dmude = getGradientShearModulus(this,material,strain,ev,ned,K)
            if strcmp(material.asympt,'vonMises')
                dmudev = 0.0;
                dmudned = -(sqrt(2.0)/2.0)*material.tauy/((material.eref+ned)^2.0);
            elseif  strcmp(material.asympt,'druckerPrager')
                dmudev = -K * material.friction / 2.0 / (material.eref + ned);
                dmudned = -material.friction * (material.sy - K * ev) / 2.0 / (material.eref + ned)^2;
            end
            % Get the derivative of the strain invariants wrt the strain
            devde = this.gradientStrainInvariantI1();
            dnedde = this.gradientNormDeviatoricStrain(strain);
            Id = this.fourthOrderSymTensor();
            % Gradient of the shear modulus (mu)
            dmude = dmudev * devde + dmudned * Id * dnedde;
        end
    end
end