%% ControlMethod_Displ Class
% This class inherits from the base class 'ControlMethod' to implement
% the displacement control method for solving nonlinear systems.
%
%% Authors
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
% * Rafael Rangel (rrangel@cimne.upc.edu)
% 
%% Class definition
classdef ControlMethod_Displ < ControlMethod
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = ControlMethod_Displ()
            this = this@ControlMethod();
        end
    end

    %% Public methods
    methods
        %------------------------------------------------------------------
        % Compute inrement of load ratio for the predicted solution
        % (first iteration of first step).
        function d_lbd0 = predictedIncrementFirst(this,~,mdl,sign,D_U,d_Up0)
            d_lbd0 = this.predictedIncrement(this,mdl,sign,1,1,0.0,0.0,D_U,d_Up0,Fref);
        end

        %------------------------------------------------------------------
        % Compute inrement of load ratio for the predicted solution
        % (first iteration).
        function d_lbd0 = predictedIncrement(~,anl,mdl,sign,J,~,~,~,~,d_Up0,~)
            d_Up0 = d_Up0(mdl.doffree);
            d_lbd0 = J * sign * anl.increment / d_Up0(anl.ctrlDof);
        end

        %------------------------------------------------------------------
        % Compute inrement of load ratio for the corrected solutions
        % (iterations to correct predicted solution).
        function d_lbd = correctedIncrement(~,anl,mdl,~,~,~,~,d_Up,d_Ur,~,~,~)
            d_Up = d_Up(mdl.doffree);
            d_Ur = d_Ur(mdl.doffree);
            d_lbd = -d_Ur(anl.ctrlDof)/d_Up(anl.ctrlDof);
        end
    end
end
