%% Anl_Transient Class
% This class inherits from the base class 'Anl' to implement the solution of
% a transient nonlinear analysis with implicit time integration schemes.
%
%% Authors
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
% 
%% Class definition
classdef Anl_Transient < Anl
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        nlscheme    = [];     % Nonlinear solution scheme
        theta       = 1.0;    % Implicit time integration scheme parameter
        ti          = 0.01;   % Initial time
        tf          = 1.0;    % Final time
        dt          = 0.001;  % Time step
        dtMax       = 0.001;  % Maximum time step
        dtMin       = 0.001;  % Minimum time step
        adaptStep   = false;  % Adaptive step size
        maxIter     = 250;    % Maximum number of iterations
        maxAttempts = 10;     % Maximum attempts to converge
        desiredIter = 5;      % Desired number of iterations
        echo        = true;   % Flag to print in the command window
    end

    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Anl_Transient(nlscheme)
            this = this@Anl('Transient');

            if strcmp(nlscheme,'Picard')
                this.nlscheme = NonlinearScheme_Picard();
            elseif strcmp(nlscheme,'Newton')
                this.nlscheme = NonlinearScheme_Newton();
            else
                disp("Error creating the Analysis object.");
                disp("Nonlinear solution scheme was not provided.");
                disp("Available options:");
                disp("   Picard");
                disp("   Newton");
                error("Error: Nonlinear solution scheme was not provided");
            end
        end
    end

    %% Public methods
    methods
        %------------------------------------------------------------------
        % Execute the transient nonlinear analysis, handle time-stepping,
        % convergence checks, and update the model state.
        % TODO: Improve it to avoid unnecessary iterations
        function run(this,mdl)
            disp("*** Performing transient nonlinear analysis...")

            % Initialize model object
            mdl.preComputations();

            % Initialize analysis parameters
            t    = this.ti + this.dt;
            t0   = this.ti;
            step = 1;

            % Initialize solution vector
            X  = mdl.U;
            dx = zeros(mdl.ndof);

            % Start transient process
            attempt = 1;
            brokenStep = false;
            while (t0 < this.tf)
                if this.echo
                    fprintf("\t Time: %12.5f s \n", t);
                end

                % Update transient solution
                XOld = X;

                % Start iterative process
                convFlg = false;
                attemptOld = attempt;
                attempt = 1;

                while attempt < this.maxAttempts
                    iter = 1;

                    while true
                        % Compute model global matrices
                        [A,b] = mdl.getLinearSystem(X,XOld,this.nlscheme,this.dt);

                        % Apply Dirichlet boundary conditions
                        [A,b] = mdl.applyDirichletBC(A,b,X,this.nlscheme);

                        % Update variables
                        [X,dx] = this.nlscheme.eval(A,b,X,dx,mdl.doffree,iter);

                        % Check convergence
                        convFlg = this.nlscheme.convergence(X,XOld,dx,b,mdl.doffree,iter,this.echo);
                        if convFlg == true
                            break;
                        end
                        
                        % Check for NaN
                        if (any(isnan(dx)) || any(isnan(X)))
                            break;
                        end

                        % Check maximum number of iterations
                        iter = iter + 1;
                        if (iter > this.maxIter)
                            break
                        end
                    end

                    % Check convergence
                    if convFlg == true
                        break;
                    end

                    % Reduce time step 
                    this.dt = max(this.dt/4.0, this.dtMin);

                    % Clean previous attempt
                    X = XOld;

                    % Update attempt counter
                    attempt = attempt + 1;
                end

                if convFlg == false
                    disp("Solution did not converge!");
                    break;
                else
                    mdl.updateStateVar();
                end

                % Update time step
                if (this.adaptStep == true) && (attempt == 1) && (brokenStep == false) && (attemptOld == 1)
                    fstep = (this.desiredIter/iter)^(0.25);
                    this.dt = max(min(fstep * this.dt, this.dtMax),this.dtMin);
                end

                % Update time
                t0 = t;
                if (t + this.dt) > this.tf
                    this.dt = this.tf - t;
                    if (abs(t - this.tf)< 1.0e-15) && (abs(this.dt) < 1.0e-12), break, end
                end
                t = t + this.dt;
                step = step + 1;
            end

            % Update state variables
            mdl.updateStateVar();

            % Save final result
            mdl.U = X;
            for i = 1:mdl.nelem
                gle = mdl.element(i).type.gle;
                mdl.element(i).type.ue = mdl.U(gle);
            end

            disp("*** Analysis completed!");
        end

        %------------------------------------------------------------------
        % Configure the transient solver with specified parameters.
        function setUpTransientSolver(this,ti,dt,tf,dtMax,dtMin,adaptStep)
            if nargin == 4
                dtMax = dt;
                dtMin = dt;
                adaptStep = false;
            end
            this.ti = ti;
            this.dt = dt;
            this.tf = tf;
            this.adaptStep = adaptStep;
            this.dtMax = dtMax;
            this.dtMin = dtMin;
        end

        %------------------------------------------------------------------
        % Enable or disable scaling of the linear system
        function setScaleLinearSystem(this,flag)
            this.nlscheme.scaleLinearSystem = flag;
        end

        %------------------------------------------------------------------
        % Enable or disable relaxation for the Picard nonlinear solution scheme.
        function setPicardRelaxation(this,flag)
            this.nlscheme.applyRelaxation = flag;
        end

        %------------------------------------------------------------------
        % Enable or disable normalization of the error for relative convergence criteria.
        function setRelativeConvergenceCriteria(this,flag)
            this.nlscheme.normalizeError = flag;
        end

        %------------------------------------------------------------------
        % Prints the solution vector 'X' for each node in the model.
        function printStep(~,X,mdl)
            for i = 1:mdl.nnodes
                fprintf("  %4d: \t", i);
                for j = 1:mdl.ndof_nd
                    fprintf("  %8.4f ", X(mdl.ID(i,j)));
                end
                fprintf("\n");
            end
        end
    end
end
