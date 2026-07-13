%% NonlinearScheme_Newton Class
% This class inherits from the base class 'NonlinearScheme' to implement
% a fully implicit time integration scheme using the Newton-Raphson method for solving nonlinear systems.
%
%% Authors
% * Danilo Cavalcanti
% 
%% Class definition
classdef NonlinearScheme_Newton < NonlinearScheme
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = NonlinearScheme_Newton()
            this = this@NonlinearScheme();
        end
    end

    %% Public methods
    methods
        %------------------------------------------------------------------
        % Assemble the Jacobian matrix and the residual vector.
        function [J,r] = assembleLinearSystem(~,C,K,fi,fe,dfidx,x,xOld,dt)
            % Residual vector
            r = fi + K * x + C * (x - xOld) / dt - fe;

            % Jacobian matrix
            J = K + dfidx + C / dt;
        end

        %------------------------------------------------------------------
        % Apply boundary conditions to the right-hand side of the system.
        function bf = applyBCtoRHS(~,~,b,~,doffree,~)
            bf = b(doffree);
        end

        %------------------------------------------------------------------
        % Add nodal forces to the right-hand side vector.
        function b = addNodalForces(~,b,fe)
            b = b - fe;
        end

        %------------------------------------------------------------------
        % Evaluate the solution increment and updates the solution vector.
        function [X,dx] = eval(this,J,r,X,~,freedof,~)
            % Compute increment of variables
            if this.scaleLinearSystem
                dx = this.solveScaledSystem(J,-r);
            else
                dx = -J\r;
            end

            % Update variables
            X(freedof) = X(freedof) + dx;
        end
        
        %------------------------------------------------------------------
        % Scale and solve the linear system A * x = b
        function x = solveScaledSystem(~,A,b)
            % Get size of the system
            n  = size(A,1);
            % Get the absolute value of the diagonal terms in A
            d  = abs(diag(A));
            % Tiny number in A's precision
            epsd = eps(class(full(1)));                
            % Scale factors
            s  = 1 ./ sqrt(max(d, epsd));
            % Scaling matrix
            S  = spdiags(s, 0, n, n);
            % Sym. scaled system
            As = S * A * S;                             
            bs = S * b;
            % Solve
            y  = As \ bs;                               
            % Unscale back
            x  = S * y;                                 
        end

        %------------------------------------------------------------------
        % Check for convergence of the nonlinear scheme.
        function convFlg = convergence(this,~,XOld,dx,r,~,iter,echo)
            normXOld = norm(XOld);
            if normXOld < 1.0e-16, normXOld = 1.0; end
            if echo
                fprintf("\t\t iter.: %3d , ||R|| = %7.3e  , ||dx||/||X0|| = %7.3e \n",iter,norm(r),norm(dx)/normXOld);
            end    
            if ((norm(r) < this.tol) || (norm(dx)/normXOld) < this.tol) && (iter > 1)
                convFlg = true;
            else
                convFlg = false;
            end
        end
    end
end
