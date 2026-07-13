%% ControlMethod_ArcSph Class
% This class inherits from the base class 'ControlMethod' to implement
% the arc length (spherical version) control method for solving nonlinear systems.
%
%% Authors
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
% * Rafael Rangel (rrangel@cimne.upc.edu)
% 
%% Class definition
classdef ControlMethod_ArcSph < ControlMethod
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = ControlMethod_ArcSph()
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
            d_lbd0 = sign * J * sqrt((D_U'*D_U + D_lbd^2*(Pref'*Pref)) / (d_Up0'*d_Up0 + Pref'*Pref));
        end

        %------------------------------------------------------------------
        % Compute inrement of load ratio for the corrected solutions
        % (iterations to correct predicted solution).
        function d_lbd = correctedIncrement(~,~,mdl,~,D_lbd,~,~,d_Up,d_Ur,D_U,Pref,~)
            d_Up = d_Up(mdl.doffree);
            d_Ur = d_Ur(mdl.doffree);
            D_U = D_U(mdl.doffree);
            Pref = Pref(mdl.doffree);
            a = d_Up'*d_Up + Pref'*Pref;
            b = d_Up'*(d_Ur + D_U) + D_lbd*(Pref'*Pref);
            c = d_Ur'*(d_Ur + 2*D_U);
            s = sign(D_U'*d_Up);
            d_lbd = -b/a + s*sqrt((b/a)^2 - c/a);
        end
    end
end
