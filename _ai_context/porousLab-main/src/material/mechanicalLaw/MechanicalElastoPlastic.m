%% MechanicalElastoPlastic Class
% This class implements the common interface for elasto-plastic
% constitutive laws. It inherits from the _MechanicalLinearElastic_ class
% and defines the stress integration data required by the implicit
% return-mapping algorithm.
%
%% Methods
% * *eval*: Computes the stress vector and the constitutive matrix.
% * *yieldCondition*: Computes the yield function value.
% * *yieldStressGradient*: Computes the gradient of the yield function 
%                          with respect to the stress vector.
% * *yieldStateGradient*: Computes the gradient of the yield function
%                         with respect to the state variables.
% * *flowVector*: Computes the flow direction vector (plastic strain 
%                 direction).
% * *flowStressGradient*: Computes the gradient of the flow vector with
%                         respect to the stress vector.
% * *flowStateGradient*: Computes the gradient of the flow vector with
%                        respect to the state variables.
% * *stateEvolution*: Computes the hardening/softening law.
% * *stateStressGradient*: Computes the gradient of the hardening/softening
%                          law with respect to the stress vector.
% * *stateStateGradient*: Computes the gradient of the hardening/softening
%                         law with respect to the state variables.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class Definition
classdef MechanicalElastoPlastic < MechanicalLinearElastic  
    properties (SetAccess = public, GetAccess = public)
        returnMappingMaxIter = 100;
        returnYieldConditionTol = 1.0e-8;
        returnNormResidualFlowRuleTol = 1.0e-8;
        returnNormResidualHardeningTol = 1.0e-8;
        nPlasticInternalVars = 0;
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MechanicalElastoPlastic()
            this = this@MechanicalLinearElastic();
        end
    end
    %% Abstract methods
    methods(Abstract)

        % Yield surface and its gradients ---------------------------------
        f     = yieldCondition(this, material, ip, stress, state);
        dfds  = yieldStressGradient(this, material, ip, stress, state);
        dfda  = yieldStateGradient(this, material, ip, stress, state);

        % Flow vector and its gradients -----------------------------------
        n     = flowVector(this, material, ip, stress, state);
        dnds  = flowStressGradient(this, material, ip, stress, state);
        dnda  = flowStateGradient(this, material, ip, stress, state);

        % Internal/state variables evolution and its gradients ------------
        h     = stateEvolution(this, material, ip, stress, state);
        dhds  = stateStressGradient(this, material, ip, stress, state);
        dhda  = stateStateGradient(this, material, ip, stress, state);

    end
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Compute the stress vector and the constitutive matrix
        function [stress,Dt] = eval(this,material,ip)
            [stress,Dt] = this.implicitReturnMapping(material,ip);
        end

        %------------------------------------------------------------------
        % Compute the stress vector and the constitutive matrix with the
        % fully implicit return mapping algorithm
        function [stress,Dt] = implicitReturnMapping(this,material,ip)

            % Constitutive matrix
            De = this.elasticConstitutiveMatrix(material,ip);
            Ce = this.elasticFlexibilityMatrix(material,ip);
            Dt = De;

            % Trial stress vector
            stress = De * (ip.strain - ip.strainOld) + ip.stressOld;

            % Previous plastic variables
            epOld    = ip.plasticstrainOld;
            stateOld = ip.statevarOld;
            state    = stateOld;
            nStress  = length(stress);
            nState   = length(stateOld);

            % Evaluate the yield condition at the trial state
            ip.plasticstrain = epOld;
            ip.statevar      = stateOld;
            f = this.yieldCondition(material,ip,stress,state);
            yieldTol = this.returnYieldConditionTol * max(1.0,norm(stress));

            % Elastic step
            if f < yieldTol
                return
            end

            % Initialize variables for the return mapping
            lambda = 0.0;
            iter = 1;
            converged = false;

            % Return mapping: closest point projection
            while iter <= this.returnMappingMaxIter

                % Yield function, flow vector and hardening/softening law
                f = this.yieldCondition(material,ip,stress,state);
                n = this.flowVector(material,ip,stress,state);
                h = this.stateEvolution(material,ip,stress,state);
                yieldTol = this.returnYieldConditionTol * max(1.0,norm(stress));

                % Residuals of the plastic flow and hardening laws
                r1 = Ce * stress - ip.strain + epOld + lambda * n;
                r2 = -state + stateOld + lambda * h;

                % Check convergence
                if (abs(f) <= yieldTol) && ...
                   (norm(r1) <= this.returnNormResidualFlowRuleTol) && ...
                   (norm(r2) <= this.returnNormResidualHardeningTol)
                    converged = true;
                    break
                end

                % Gradients of the yield function
                dfds = this.yieldStressGradient(material,ip,stress,state);
                dfda = this.yieldStateGradient(material,ip,stress,state);

                % Gradients of the flow rule vector
                dnds = this.flowStressGradient(material,ip,stress,state);
                dnda = this.flowStateGradient(material,ip,stress,state);

                % Gradients of the hardening/softening law
                dhds = this.stateStressGradient(material,ip,stress,state);
                dhda = this.stateStateGradient(material,ip,stress,state);

                % Auxiliary system from the first two residual blocks
                K = [Ce + lambda * dnds, lambda * dnda;
                     lambda * dhds, -eye(nState) + lambda * dhda];
                rho = [r1; r2];
                e = [n; h];
                g = [dfds; dfda];

                % Increment of the plastic multiplier
                dlambda = (f - g' * (K \ rho)) / (g' * (K \ e));

                % Stress and internal-variable corrections
                delta = -K \ (rho + dlambda * e);
                stress = stress + delta(1:nStress);
                state  = state  + delta(nStress+1:end);
                lambda = lambda + dlambda;

                % Update iteration counter
                iter = iter + 1;

            end

            % Check convergence of the local Newton iteration
            if ~converged
                error(['Error: implicit return mapping did not converge. ', ...
                       'f = %.6e, yieldTol = %.6e, norm(r1) = %.6e, norm(r2) = %.6e, lambda = %.6e, iter = %d'], ...
                       f,yieldTol,norm(r1),norm(r2),lambda,iter);
            end

            % Store the accepted state variables
            ip.statevar = state;
            ip.plasticstrain = ip.strain - Ce * stress;

            % Compute algorithmic tangent constitutive tensor
            n = this.flowVector(material,ip,stress,state);
            h = this.stateEvolution(material,ip,stress,state);
            dfds = this.yieldStressGradient(material,ip,stress,state);
            dfda = this.yieldStateGradient(material,ip,stress,state);
            dnds = this.flowStressGradient(material,ip,stress,state);
            dnda = this.flowStateGradient(material,ip,stress,state);
            dhds = this.stateStressGradient(material,ip,stress,state);
            dhda = this.stateStateGradient(material,ip,stress,state);

            K = [Ce + lambda * dnds, lambda * dnda;
                 lambda * dhds, -eye(nState) + lambda * dhda];
            e = [n; h];
            g = [dfds; dfda];
            rhs = [eye(nStress); zeros(nState,nStress)];

            Krhs = K \ rhs;
            Ke   = K \ e;
            dsol = Krhs - Ke * ((g' * Krhs) / (g' * Ke));
            Dt = dsol(1:nStress,:);

        end
    end
    methods (Static)
        %------------------------------------------------------------------
        % Flag to impose if the material is elasto-plastic or not
        function flag = isElastoPlastic()
            flag = true;
        end
    end
end
