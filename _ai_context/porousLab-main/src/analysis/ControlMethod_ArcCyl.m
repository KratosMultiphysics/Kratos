%% ControlMethod_ArcCyl Class
% This class inherits from the base class 'ControlMethod' to implement
% the arc length (cylindrical version) control method for solving nonlinear systems.
%
%% Authors
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
% * Rafael Rangel (rrangel@cimne.upc.edu)
% 
%% Class definition
classdef ControlMethod_ArcCyl < ControlMethod
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = ControlMethod_ArcCyl()
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
        function d_lbd = correctedIncrement(~,~,mdl,~,~,~,~,d_Up,d_Ur,D_U,~,~)
            d_Up = d_Up(mdl.doffree);
            d_Ur = d_Ur(mdl.doffree);
            D_U = D_U(mdl.doffree);
            a = d_Up'*d_Up;
            b = d_Up'*(d_Ur + D_U);
            c = d_Ur'*(d_Ur + 2*D_U);
            s = sign(D_U'*d_Up);
            d_lbd = -b/a + s*sqrt((b/a)^2 - c/a);
        end
    end
end
