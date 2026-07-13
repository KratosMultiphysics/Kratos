%% ControlMethod_OrtResidual Class
% This class inherits from the base class 'ControlMethod' to implement
% the orthogonal residual control method for solving nonlinear systems.
%
%% Authors
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
% * Rafael Rangel (rrangel@cimne.upc.edu)
% 
%% Class definition
classdef ControlMethod_OrtResidual < ControlMethod
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = ControlMethod_OrtResidual()
            this = this@ControlMethod();
        end
    end

    %% Public methods
    methods
        %------------------------------------------------------------------
        % Compute inrement of load ratio for the predicted solution
        % (first iteration of first step).
        function d_lbd0 = predictedIncrementFirst(~,anl,~,~,~,~)
            d_lbd0 = anl.increment;
        end

        %------------------------------------------------------------------
        % Compute inrement of load ratio for the predicted solution
        % (first iteration).
        function d_lbd0 = predictedIncrement(~,~,mdl,sign,J,~,~,~,D_U,d_Up0,~)
            D_U = D_U(mdl.doffree);
            d_Up0 = d_Up0(mdl.doffree);
            d_lbd0 = sign * J * sqrt((D_U'*D_U)/(d_Up0'*d_Up0));
        end

        %------------------------------------------------------------------
        % Compute inrement of load ratio for the corrected solutions
        % (iterations to correct predicted solution).
        function d_lbd = correctedIncrement(~,~,mdl,~,~,~,~,~,~,D_U,Pref,R)
            D_U = D_U(mdl.doffree);
            Pref = Pref(mdl.doffree);
            R = R(mdl.doffree);
            d_lbd = -(R'*D_U)/(Pref'*D_U);
        end
    end
end
