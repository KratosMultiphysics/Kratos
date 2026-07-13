%% Anl_Linear Class
% This class inherits from the base class 'Anl' to implement the solution a static linear analysis.
%
%% Authors
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
% 
%% Class definition
classdef Anl_Linear < Anl
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
    end

    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Anl_Linear()
            this = this@Anl('Linear');
        end
    end

    %% Public methods
    methods
        %------------------------------------------------------------------
        % Execute the linear analysis for the given model object 'mdl'.
        function run(~,mdl)
            disp("*** Performing linear analysis...")

            % Initialize model object
            mdl.preComputations();

            % Compute global stiffness matrix
            [K,~,Fi,Fext,dFidX] = mdl.globalMatrices(mdl.U);

            % Set linear system
            if ((nnz(K)>0) && (nnz(dFidX) == 0))
                A = K;
            elseif ((nnz(K) == 0) && (nnz(dFidX) > 0))
                A = dFidX;
            elseif ((nnz(K) == 0) && (nnz(dFidX) == 0))
                error("Physics did not fill either K or dFidX.")
            end
            b = Fext(mdl.doffree) - Fi(mdl.doffree) - A(mdl.doffree,mdl.doffixed) * mdl.U(mdl.doffixed);

            % Solve linear system
            mdl.U(mdl.doffree) = A(mdl.doffree,mdl.doffree)\b;

            % Save final result
            for i = 1:mdl.nelem
                gle = mdl.element(i).type.gle;
                mdl.element(i).type.ue = mdl.U(gle);
            end

            % Call it again to update state variables
            mdl.globalMatrices(mdl.U);
            mdl.updateStateVar();

            disp("*** Analysis completed!");
        end
    end
end
