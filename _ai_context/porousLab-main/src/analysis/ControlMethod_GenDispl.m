%% ControlMethod_GenDispl Class
% This class inherits from the base class 'ControlMethod' to implement
% the generalized displacement control method for solving nonlinear systems.
%
%% Authors
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
% * Rafael Rangel (rrangel@cimne.upc.edu)
% 
%% Class definition
classdef ControlMethod_GenDispl < ControlMethod
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = ControlMethod_GenDispl()
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
        function d_lbd0 = predictedIncrement(~,anl,~,sign,J,GSP,~,~,~,~,~)
            d_lbd0 = sign * J * sqrt(abs(GSP)) * anl.increment;
        end

        %------------------------------------------------------------------
        % Compute inrement of load ratio for the corrected solutions
        % (iterations to correct predicted solution).
        function d_lbd = correctedIncrement(~,~,mdl,~,~,d_Up0,~,d_Up,d_Ur,~,~,~)
            % Extract free DOF components
            d_Up0 = d_Up0(mdl.doffree);
            d_Up  = d_Up(mdl.doffree);
            d_Ur  = d_Ur(mdl.doffree);
            d_lbd = -(d_Up0'*d_Ur)/(d_Up0'*d_Up);
        end
    end
end
