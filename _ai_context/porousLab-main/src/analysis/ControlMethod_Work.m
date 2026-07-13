%% ControlMethod_Work Class
% This class inherits from the base class 'ControlMethod' to implement
% the work control method for solving nonlinear systems.
%
%% Authors
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
% * Rafael Rangel (rrangel@cimne.upc.edu)
% 
%% Class definition
classdef ControlMethod_Work < ControlMethod
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = ControlMethod_Work()
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
        function d_lbd0 = predictedIncrement(~,~,mdl,sign,J,~,D_lbd,~,D_U,d_Up0,Pref)
            Pref = Pref(mdl.doffree);
            D_U = D_U(mdl.doffree);
            d_Up0 = d_Up0(mdl.doffree);
            d_lbd0 = sign * J * sqrt(abs((D_lbd*Pref'*D_U)/(Pref'*d_Up0)));
        end

        %------------------------------------------------------------------
        % Compute inrement of load ratio for the corrected solutions
        % (iterations to correct predicted solution).
        function d_lbd = correctedIncrement(~,~,mdl,~,~,~,~,d_Up,d_Ur,~,Pref,~)
            d_Up = d_Up(mdl.doffree);
            d_Ur = d_Ur(mdl.doffree);
            Pref = Pref(mdl.doffree);
            d_lbd = -(Pref'*d_Ur)/(Pref'*d_Up);
        end
    end
end
