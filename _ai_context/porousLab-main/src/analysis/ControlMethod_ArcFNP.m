%% ControlMethod_ArcFNP Class
% This class inherits from the base class 'ControlMethod' to implement
% the arc length (fixed normal plane version) control method for solving nonlinear systems.
%
%% Authors
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
% * Rafael Rangel (rrangel@cimne.upc.edu)
% 
%% Class definition
classdef ControlMethod_ArcFNP < ControlMethod
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = ControlMethod_ArcFNP()
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
            % Extract free DOF components
            D_U = D_U(mdl.doffree);
            d_Up0 = d_Up0(mdl.doffree);
            d_lbd0 = sign * J * sqrt((D_U'*D_U)/(d_Up0'*d_Up0));
        end

        %------------------------------------------------------------------
        % Compute inrement of load ratio for the corrected solutions
        % (iterations to correct predicted solution).
        function d_lbd = correctedIncrement(~,~,mdl,d_lbd0,~,~,d_U0,d_Up,d_Ur,~,Pref,~)
            d_U0 = d_U0(mdl.doffree);
            d_Up = d_Up(mdl.doffree);
            d_Ur = d_Ur(mdl.doffree);
            Pref = Pref(mdl.doffree);            
            d_lbd = -(d_Ur'*d_U0)/(d_Up'*d_U0 + d_lbd0*(Pref'*Pref));
        end
    end
end
