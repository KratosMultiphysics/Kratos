%% MechanicalElastoPlasticDruckerPrager Class
% This class implements the Drucker-Prager criteria for the elasto-plastic
% material law. It provides methods for evaluating stress, constitutive
% matrices, yield conditions, flow vectors, and their gradients, as well
% as handling plastic strain updates.
%
%% Methods
% * *eval*: Computes the stress vector and the constitutive matrix for the 
%           material at a given integration point. Handles both elastic 
%           and plastic steps.
% * *alternativeStressIntegration*: Implements an alternative stress 
%                                   integration algorithm for the material.
% * *getMohrCoulombCorrespondence*: Computes the Mohr-Coulomb 
%                                   correspondence parameters for the 
%                                   material.
% * *yieldCondition*: Defines the yield function based on the Drucker-
%                     Prager criteria.
% * *yieldStressGradient*: Computes the gradient of the yield function 
%                          with respect to the stress vector.
% * *flowVector*: Computes the flow vector for the plastic potential.
% * *flowStressGradient*: Computes the gradient of the flow vector with
%                         respect to the stress vector.
% * *pseudoInv*: Computes the pseudoinverse of a given matrix using SVD.
% * *stateEvolution*: Computes the hardening/softening law.
% * *stateStressGradient*: Computes the gradient of the hardening/softening
%                          law with respect to the stress vector.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class Definition
classdef MechanicalElastoPlasticDruckerPrager < MechanicalElastoPlastic  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MechanicalElastoPlasticDruckerPrager()
            this = this@MechanicalElastoPlastic();
            this.nstVar = 0;   % Hardening + Kinematic hardening
        end
    end

    %% Public methods
    methods
        %------------------------------------------------------------------
        % Compute the stress vector and the constitutive matrix
        function [stress,Dt] = eval(this,material,ip)
            if strcmp(material.stressIntAlgorithm,'implicit')
                 [stress,Dt] = eval@MechanicalElastoPlastic(this,material,ip);
            elseif strcmp(material.stressIntAlgorithm,'alternative')
                [stress,Dt] = this.alternativeStressIntegration(material,ip);
            else
                disp('Error: the given stress integration algorithm is not available');
                disp('Tags of the methods available: ''implicit'', ''alternative''');
                error('Error: stressIntAlgorithm is not available');
            end
        end

        %------------------------------------------------------------------
        % Compute the stress vector and the constitutive matrix
        function [stress,Dt] = alternativeStressIntegration(this,material,ip)

            % Constitutive matrix
            De = this.elasticConstitutiveMatrix(material,ip);
            Ce = this.elasticFlexibilityMatrix(material,ip);
            Dt = De;

            % Trial stress vector
            stress = De * (ip.strain - ip.strainOld) + ip.stressOld;

            % Evaluate the yield condition
            f = this.yieldCondition(material,ip,stress);

            % Elastic step
            if f < this.returnYieldConditionTol, return, end

            % Material parameters
            [eta, xi, etaB] = this.getMohrCoulombCorrespondence(material);
            coh    = material.cohesion;
            Id     = this.gradientI1(stress);

            % Elastic properties
            K = this.bulkModulus(material);
            G = this.shearModulus(material);

            % Stress invariants
            p = this.hydrostaticStress(stress);
            J2 = this.stressInvariantJ2(stress);

            % Deviatoric stresses
            s = this.deviatoricStress(stress);

            % Plastic multiplier
            lambda = f / (G + K * eta * etaB);

            % Stress update to the smooth part of the cone 
            factor = 1.0 - G * lambda / sqrt(J2);
            if J2 > 0.0    
                s = factor * s;
                p = p - lambda * etaB * K;
                stress = s + p * Id;
                df = this.yieldStressGradient(material,ip,stress);
                n  = this.flowVector(material,ip,stress);
                dn = this.flowStressGradient(material,ip,stress,[]);
                Psi = this.pseudoInv(Ce + lambda * dn);
                Dt  = Psi - ((Psi * n) * df' * Psi)/((df' * Psi) * n);
            else 
                % Return to the apex of the surface
                stress = (xi * coh / eta) * Id;
                Dt = zeros(4,4);
            end

            % Update the plastic strain
            ip.plasticstrain = ip.strain - Ce * stress;
        end

        %------------------------------------------------------------------
        function [eta, xi, etab] = getMohrCoulombCorrespondence(~,material)
            % Mohr-Coulomb parameters
            phi = material.frictionAngle;
            psi = material.dilationAngle;
            if strcmp(material.MCmatch,'inner')
                xi   = 6.0 * cos(phi) / (sqrt(3.0 * (3.0 + sin(phi))));
                eta  = 6.0 * sin(phi) / (sqrt(3.0 * (3.0 + sin(phi))));
                etab = 6.0 * sin(psi) / (sqrt(3.0 * (3.0 + sin(psi))));
            elseif  strcmp(material.MCmatch,'outer')
                xi   = 6.0 * cos(phi) / (sqrt(3.0 * (3.0 - sin(phi))));
                eta  = 6.0 * sin(phi) / (sqrt(3.0 * (3.0 - sin(phi))));
                etab = 6.0 * sin(psi) / (sqrt(3.0 * (3.0 - sin(psi))));
            elseif strcmp(material.MCmatch,'planestrain')
                xi   = 3.0 / (sqrt(9.0 + 12.0*tan(phi)*tan(phi)));
                eta  = 3.0 * tan(phi) / (sqrt(9.0 + 12.0*tan(phi)*tan(phi)));
                etab = 3.0 * tan(psi) / (sqrt(9.0 + 12.0*tan(psi)*tan(psi)));
            else
                xi   = 6.0 * cos(phi) / (sqrt(3.0 * (3.0 + sin(phi))));
                eta  = 6.0 * sin(phi) / (sqrt(3.0 * (3.0 + sin(phi))));
                etab = 6.0 * sin(psi) / (sqrt(3.0 * (3.0 + sin(psi))));
            end
        end

        %------------------------------------------------------------------
        % Yield function definition
        function f = yieldCondition(this,material,~,stress,~)
            % Material parameters
            [eta, xi] = this.getMohrCoulombCorrespondence(material);
            % Stress invariants
            p  = this.hydrostaticStress(stress);
            J2 = this.stressInvariantJ2(stress);
            % Yield surface
            f = sqrt(J2) + eta * p - xi * material.cohesion;
        end

        %------------------------------------------------------------------
        % Gradient of the yield function wrt to the stress vector
        function df = yieldStressGradient(this,material,ip,stress,~)
            % Material parameters
            eta = this.getMohrCoulombCorrespondence(material);
            % Deviatoric stress invariant
            J2 = this.stressInvariantJ2(stress);
            % Stress invariants gradients
            dI1 = this.gradientI1(stress);
            dJ2 = this.gradientJ2(stress);
            % Derivatives of the yield surface wrt to the invariants
            dfdI1 = eta / 3.0;
            if J2 > 0.0
                dfdJ2 = 1.0 /(2.0 * sqrt(J2));
            else
                dfdJ2 = 0.0;
            end
            % Yield surface gradient
            df = dfdI1 * dI1 + dfdJ2 * dJ2;
            if strcmp(ip.anm,'PlaneStress')
                df(3) = 0.0;
            end
        end

        %------------------------------------------------------------------
        % Gradient of the yield function wrt to the state variables vector
        function dfda = yieldStateGradient(~,~,~,~,~)
            dfda = zeros(0,1);
        end

        %------------------------------------------------------------------
        % Flow vector
        function n = flowVector(this,material,ip,stress,~)
            % Material parameters
            [~,~,etaB] = this.getMohrCoulombCorrespondence(material);
            % Deviatoric stress invariant
            J2 = this.stressInvariantJ2(stress);
            % Stress invariants gradients
            dI1 = this.gradientI1(stress);
            dJ2 = this.gradientJ2(stress);
            % Derivatives of the yield surface wrt to the invariants
            dfdI1 = etaB / 3.0;
            if J2 > 0.0
                dfdJ2 = 1.0 /(2.0 * sqrt(J2));
            else
                dfdJ2 = 0.0;
            end
            % Yield surface gradient
            n = dfdI1 * dI1 + dfdJ2 * dJ2;
            if strcmp(ip.anm,'PlaneStress')
                n(3) = 0.0;
            end
        end

        %------------------------------------------------------------------
        % Flow vector gradient
        function dn = flowStressGradient(this,~,ip,stress,~)
            % Deviatoric stress invariant 
            J2   = this.stressInvariantJ2(stress);
            dJ2  = this.gradientJ2(stress);
            d2J2 = this.hessianJ2();
            % Derivatives of the yield surface wrt to the invariants
            if J2 > 0.0
                dfdJ2 = 0.5 * sqrt(1.0 / J2);
                d2fdJ2 = -0.25 * sqrt(1.0 / J2 / J2 / J2);
            else
                dfdJ2 = 0.0;
                d2fdJ2 = 0.0;
            end
            % Yield surface gradient
            dn = d2fdJ2 * (dJ2 * dJ2') + dfdJ2 * d2J2;
            if strcmp(ip.anm,'PlaneStress')
                dn(3,:) = 0.0;
                dn(:,3) = 0.0;
                dn(3,3) = 1.0;
            end
        end

        %------------------------------------------------------------------
        % Flow vector gradient wrt to the state variables vector
        function dnda = flowStateGradient(~,~,ip,~,~)
            dnda = zeros(ip.nVar,0);
        end

        %------------------------------------------------------------------
        % Computes the pseudoinverse of a given matrix using SVD
        function Ai = pseudoInv(~,A)
            % Assume A is your input matrix
            [U, S, V] = svd(A);  % Compute the SVD of A
            
            % Invert the non-zero singular values in S
            S_inv = zeros(size(S'));  % Initialize a matrix for the pseudoinverse of S
            tolerance = 1e-10;  % A tolerance for small singular values
            for i = 1:min(size(S))
                if S(i, i) > tolerance
                    S_inv(i, i) = 1 / S(i, i);  % Inverse of the singular value
                end
            end
            
            % Compute the pseudoinverse of A
            Ai = V * S_inv * U';
        end

        %------------------------------------------------------------------
        % Internal/state variables evolution
        function h = stateEvolution(~,~,~,~,~)
            h = zeros(0,1);
        end

        %------------------------------------------------------------------
        % Gradient of the internal/state variables law wrt to the stress vector
        function dhds = stateStressGradient(~,~,ip,~,~)
            dhds = zeros(0,ip.nVar);
        end

        %------------------------------------------------------------------
        % Gradient of the internal/state variables law wrt to the state variables
        function dhda = stateStateGradient(~,~,~,~,~)
            dhda = zeros(0,0);
        end

    end
end
