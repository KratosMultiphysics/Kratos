%% MechanicalElastoPlasticVonMises Class
% This class implements an elasto-plastic constitutive law based on the 
% Von Mises yield criterion. It extends the _MechanicalElastoPlastic_ 
% base class and provides methods for defining the yield condition, 
% flow rule, and hardening behavior.
%
%% Methods
% * *yieldCondition*: Computes the yield condition based on the Von Mises 
%                     stress and the yield stress.
% * *yieldStressGradient*: Computes the gradient of the yield function 
%                          with respect to the stress vector.
% * *flowVector*: Computes the flow direction vector based on the 
%                 deviatoric stress.
% * *flowStressGradient*: Computes the gradient of the flow vector with
%                         respect to the stress vector.
% * *stateEvolution*: Returns the scalar isotropic hardening law.
% * *stateStressGradient*: Returns the gradient of the hardening law with
%                          respect to the stress vector.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class Definition
classdef MechanicalElastoPlasticVonMises < MechanicalElastoPlastic  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MechanicalElastoPlasticVonMises()
            this = this@MechanicalElastoPlastic();
            this.nstVar = 1;   % Hardening
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        % Yield function definition
        function f = yieldCondition(this,material,ip,stress,state)
            if nargin < 5
                state = ip.statevar;
            end
            sVM = this.vonMisesStress(stress);
            sy = material.sy0 + material.Kp * state(1);
            f  = sVM - sy;
        end

        %------------------------------------------------------------------
        % Gradient of the yield function wrt to the stress vector
        function dfds = yieldStressGradient(this,~,~,stress,~)
            dfds = this.vonMisesStressGradient(stress);
        end

        %------------------------------------------------------------------
        % Gradient of the yield function wrt to the state variables vector
        function dfda = yieldStateGradient(~,material,~,~,~)
            dfda = -material.Kp;
        end

        %------------------------------------------------------------------
        % Flow vector
        function n = flowVector(this,material,ip,stress,state)
            n = this.yieldStressGradient(material,ip,stress,state);
        end

        %------------------------------------------------------------------
        % Flow vector gradient wrt to the stress vector
        function dnds = flowStressGradient(this,~,~,stress,~)
            dnds = this.vonMisesStressHessian(stress);
        end

        %------------------------------------------------------------------
        % Flow vector gradient wrt to the state variables vector
        function dnda = flowStateGradient(~,~,ip,~,~)
            dnda = zeros(ip.nVar,1);
        end

        %------------------------------------------------------------------
        % Internal/state variables evolution
        function h = stateEvolution(~,~,~,~,~)
            h = 1.0;
        end

        %------------------------------------------------------------------
        % Gradient of the internal/state variables law wrt to the stress vector
        function dhds = stateStressGradient(~,~,ip,~,~)
            dhds = zeros(1,ip.nVar);
        end

        %------------------------------------------------------------------
        % Gradient of the internal/state variables law wrt to the state variables
        function dhda = stateStateGradient(~,~,~,~,~)
            dhda = 0.0;
        end

    end
end
