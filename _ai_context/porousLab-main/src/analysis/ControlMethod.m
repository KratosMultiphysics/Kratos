%% ControlMethod Class
% This is an abstract class that defines a control method object.
% for nonlinear quasi-static analysis.
%
%% Authors
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
% * Rafael Rangel (rrangel@cimne.upc.edu)
%
%% Class definition
classdef ControlMethod < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
    end

    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = ControlMethod()
            return;
        end
    end

    %% Abstract methods
    methods (Abstract)
        %------------------------------------------------------------------
        % Compute inrement of load ratio for the predicted solution
        % (first iteration of first step).
        d_lbd0 = predictedIncrementFirst(this,anl,mdl,sign,D_U,d_Up0);

        %------------------------------------------------------------------
        % Compute inrement of load ratio for the predicted solution
        % (first iteration).
        d_lbd0 = predictedIncrement(~,anl,mdl,sign,J,GSP,D_lbd,d_lbd0,D_U,d_Up0,Pref);

        %------------------------------------------------------------------
        % Compute inrement of load ratio for the corrected solutions
        % (iterations to correct predicted solution).
        d_lbd = correctedIncrement(~,anl,mdl,d_lbd0,D_lbd,d_Up0,d_U0,d_Up,d_Ur,D_U,Pref,R);
    end
end
