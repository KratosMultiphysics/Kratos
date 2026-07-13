%% MechanicalLinearElastic Class
% This class defines an linear elastic stress-strain constitutive law. The 
% class provides methods to compute the stress vector, constitutive 
% matrix, elastic constants (shear modulus and bulk modulus), and elastic
% tensors (constitutive and flexibility matrices). It also includes a 
% method to determine if the material is elasto-plastic.
%
%% Methods
% * *eval*: Computes the stress vector and the constitutive matrix based 
%           on the material properties and integration point data.
% * *shearModulus*: Computes the shear modulus using Young's modulus and 
%                   Poisson's ratio.
% * *bulkModulus*: Computes the bulk modulus using Young's modulus and 
%                  Poisson's ratio.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class Definition
classdef MechanicalLinearElastic < MechanicalLaw  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MechanicalLinearElastic()
            this = this@MechanicalLaw();
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        % Compute the stress vector and the constitutive matrix
        function [stress,De] = eval(this,material,ip)

            % Constitutive matrix
            De = this.elasticConstitutiveMatrix(material,ip);

            % Stress vector
            stress = De * (ip.strain - ip.strainOld) + ip.stressOld;

        end
        %% Elastic constants

        %------------------------------------------------------------------
        % Computes the shear modulus of the material
        function G = shearModulus(~,material)
            G = material.Young / (2.0 * (1.0 + material.nu));
        end

        %------------------------------------------------------------------
        % Computes the bulk modulus of the material
        function K = bulkModulus(~,material)
            K = material.Young / (3.0 * (1.0 - 2.0*material.nu));
        end 
    end

    %% Public methods
    methods (Static)
        function flag = isElastoPlastic()
            flag = false;
        end
        %% Elastic tensors
        
        %------------------------------------------------------------------
        % Compute the elastic constitutive matrix
        function De = elasticConstitutiveMatrix(material,ip)

            % Elastic material properties
            E  = material.Young;
            nu = material.nu;

            if strcmp(ip.anm,'PlaneStress')

                c = E/(1.0 - (nu*nu));
                De = [  c   ,   c*nu , 0.0  ,  0.0;
                       c*nu ,    c   , 0.0  ,  0.0;
                       0.0  ,   0.0  , 1.0  ,  0.0;
                       0.0  ,   0.0  , 0.0  , c*(1-nu)/2.0 ];

            else

                De = [ 1.0-nu ,   nu   ,   nu   ,    0.0;
                         nu   , 1.0-nu ,   nu   ,    0.0;
                         nu   ,  nu    , 1.0-nu ,    0.0;
                        0.0   ,  0.0   ,   0.0  , (1-2.0*nu)/2.0 ];

                De = De * E/(1.0 + nu)/(1.0 - 2.0*nu);

            end
        end

        %------------------------------------------------------------------
        % Compute the elastic flexibility matrix
        function Ce = elasticFlexibilityMatrix(material,ip)

            % Elastic material properties
            E  = material.Young;
            nu = material.nu;

            if strcmp(ip.anm,'PlaneStress')

                Ce = [  1.0/E ,  -nu/E  ,  0.0  ,  0.0;
                        -nu/E ,  1.0/E  ,  0.0  ,  0.0;
                         0.0  ,  0.0    ,  1.0  ,  0.0;
                         0.0  ,  0.0    ,  0.0  , 2.0*(1+nu)/E ];

            elseif strcmp(ip.anm,'PlaneStrain')

                Ce = [  1.0/E ,  -nu/E  ,  -nu/E ,  0.0;
                        -nu/E ,  1.0/E  ,  -nu/E ,  0.0;
                        -nu/E ,  -nu/E  ,  1.0/E ,  0.0;
                         0.0  ,  0.0    ,  0.0   , 2.0*(1+nu)/E ];
            else
                Ce = [];
            end
        end   
    end
end