%% MechanicalCohesiveLinearElastic Class
% This class implements a linear elastic cohesive law for mechanical 
% materials. It provides methods to compute the stress vector and the 
% constitutive matrix based on the material properties and the strain 
% state at integration points.
%
%% Methods
% * *eval*: Computes the stress vector and the constitutive matrix for the 
%           given material and integration point.
% * *isElastoPlastic*: Static method that indicates that the material is 
%                      not elasto-plastic.
% * *elasticConstitutiveMatrix*: Static method that computes the elastic 
%                                constitutive matrix based on the material 
%                                properties and the strain state at the 
%                                integration point.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class Definition
classdef MechanicalCohesiveLinearElastic < handle  
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        nstVar = 1;   % Number of state variables
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MechanicalCohesiveLinearElastic()
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

        %------------------------------------------------------------------
        % Compute the elastic constitutive matrix
        function De = elasticConstitutiveMatrix(this,material,ip)

            % Elastic material properties
            ks = material.shearStiffness;
            kn = this.closureModel(material,ip);
            
            % Assemble constitutive matrix
            De = [ ks   0.0;
                   0.0   kn ];
        end 
    end
    %% Public methods
    methods (Static)
        %------------------------------------------------------------------
        % Flag to return that the material is not elasto-plastic
        function flag = isElastoPlastic()
            flag = false;
        end
        %------------------------------------------------------------------
        function initializeAperture(material, ip)
            ip.statevar(1) = material.initialAperture;
            ip.statevarOld(1) = material.initialAperture;
        end
        
        %------------------------------------------------------------------
        % Compute the normal stiffness based on the closure model
        function kn = closureModel(material,ip)
            w0 = material.initialAperture;
            kn = material.normalStiffness;
            wn = ip.strain(2) + w0;
            if wn < 0
                if strcmp(material.contactPenalization,'constant')
                    kn = kn * 1.0E4;
                elseif strcmp(material.contactPenalization,'bartonBandis')
                    wmax = material.maximumClosure;
                    wn = max(wn,-wmax + 1.0E-6);
                    kn = wmax * wmax * kn / ((wn + wmax) * (wn + wmax));
                end
            end

        end
    end
end