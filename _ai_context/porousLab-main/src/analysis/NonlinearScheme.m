%% NonlinearScheme Class
% This is an abstract class that defines a nonlinear scheme object.
% It serves as a base class for implementing nonlinear solution schemes in numerical analysis.
% The class provides a set of abstract methods that must be implemented by subclasses,
% as well as some public properties and methods for managing convergence tolerance and error normalization.
%
%% Authors
% * Danilo Cavalcanti
% 
%% Class definition
classdef NonlinearScheme < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        tol = 1.0e-5;
        normalizeError = false;
        scaleLinearSystem = false;
    end

    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = NonlinearScheme()
            return;
        end
    end

    %% Abstract methods
    methods (Abstract)
        %------------------------------------------------------------------
        % Assemble the linear system for the nonlinear problem.
        [A,b] = assembleLinearSystem(this,C,K,fi,fe,dfidx,x,xOld,dt);

        %------------------------------------------------------------------
        % Apply boundary conditions to the right-hand side of the system.
        bf = applyBCtoRHS(this,A,b,x,doffree,doffixed);

        %------------------------------------------------------------------
        % Add nodal forces to the right-hand side vector.
        b = addNodalForces(this,b,fe);

        %------------------------------------------------------------------
        % Evaluates the solution increment and updates the solution vector.
        [X,dx] = eval(this,J,r,X,dx,freedof,iter);

        %------------------------------------------------------------------
        % Check for convergence of the nonlinear scheme.
        convFlg = convergence(this,X,XOld,dx,b,doffree,iter);
    end

    %% Public methods
    methods
        %------------------------------------------------------------------
        % Set the convergence tolerance.
        function setConvergenceTolerance(this,tol)
            this.tol = tol;
        end
    end
end
