%% Anl_NonlinearQuasiStatic Class
% This class inherits from the base class 'Anl' to implement the solution of
% a quasi-static nonlinear incremental-iterative analysis using different control methods.
%
%% Authors
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
% * Rafael Rangel (rrangel@cimne.upc.edu)
% 
%% Class definition
classdef Anl_NonlinearQuasiStatic < Anl
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        method        = [];      % Solution method
        adjustStep    = false;   % Flag for type of increment size adjustment
        increment     = 0.1;     % Initial increment of load ratio
        max_increment = 0.5;     % Maximum increment of load ratio
        max_lratio    = 1.0;     % Limit value of load ratio
        max_step      = 10;      % Maximum number of steps
        max_iter      = 10;      % Maximum number of iterations in each step
        trg_iter      = 3;       % Desired number of iterations in each step
        tol           = 1.0e-5;  % Numerical tolerance for convergence
        ctrlDof       = 1;       % Control DOF (for displacement control method)
        plotNd        = 1;       % Node whose DOF will be plotted
        plotDof       = 1;       % Node's DOF (ux,uy) that will plotted
        Uplot         = [];      % Matrix of nodal displacement vectors of all steps/modes
        lbdplot       = [];      % Vector of load ratios of all steps
        echo          = true;    % Flag to print in the command window
    end

    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Anl_NonlinearQuasiStatic(method)
            this = this@Anl('NonlinearQuasiStatic');

            if strcmp(method,'LoadControl')
                this.method = ControlMethod_Load();
            elseif strcmp(method,'DisplacementControl')
                this.method = ControlMethod_Displ();
            elseif strcmp(method,'WorkControl')
                this.method = ControlMethod_Work();
            elseif strcmp(method,'ArcLengthFNPControl')
                this.method = ControlMethod_ArcFNP();
            elseif strcmp(method,'ArcLengthUNPControl')
                this.method = ControlMethod_ArcUNP();
            elseif strcmp(method,'ArcLengthCylControl')
                this.method = ControlMethod_ArcCyl();
            elseif strcmp(method,'ArcLengthSPHControl')
                this.method = ControlMethod_ArcSph();
            elseif strcmp(method,'MinimumNorm')
                this.method = ControlMethod_MinNorm();
            elseif strcmp(method,'OrthogonalResidual')
                this.method = ControlMethod_OrtResidual();
            elseif strcmp(method,'GeneralizedDisplacement')
                this.method = ControlMethod_GenDispl();
            else
                disp("Error creating the Analysis object.");
                disp("Control method was not provided.");
                disp("Available options:");
                disp("   LoadControl");
                disp("   ControlMethod_Displ");
                disp("   WorkControl");
                disp("   ArcLengthFNPControl");
                disp("   ArcLengthUNPControl");
                disp("   ArcLengthCylControl");
                disp("   ArcLengthSPHControl");
                disp("   MinimumNorm");
                disp("   OrthogonalResidual");
                disp("   GeneralizedDisplacement");
                error("Error: Control method was not provided");
            end
        end
    end

    %% Public methods
    methods
        %------------------------------------------------------------------
        % Execute the nonlinear analysis process, handle iterative
        % steps, convergence checks, and state variables updates.
        function run(this,mdl)
            disp("*** Performing quasi-static nonlinear analysis...")

            % Initialize model object
            mdl.preComputations();

            % Initialize results
            this.lbdplot = zeros(this.max_step+1,1);
            this.Uplot   = zeros(this.max_step+1,1);

            % Initialize data for first step
            step = 0;  % step number
            lbd  = 0;  % total load ratio (lambda)
            sign = 1;  % sign of predicted increment of load ratio

            % Initialize vector of total nodal displacements
            U = mdl.U;

            % Initialize vector of total increment displacement
            D_U = zeros(mdl.ndof,1);

            % Start incremental process
            while (step < this.max_step)
                step = step + 1;
                if this.echo
                    fprintf("\n Step: %-4d \n", step);
                end

                % Tangent stiffness matrix
                [~,~,~,Fref,K] = mdl.globalMatrices(U);

                % Tangent increment of displacements for predicted solution
                d_Up0 = this.solveSystem(mdl,K,Fref,U);

                if (step == 1)
                    % Initial increment of load ratio for predicted solution
                    d_lbd0 = this.method.predictedIncrementFirst(this,mdl,sign,D_U,d_Up0);

                    % Set previous tangent increment of displacements as current increment
                    d_Up0_old = d_Up0;

                    % Store squared value of the norm of tangent increment of displacements
                    n2 = d_Up0(mdl.doffree)' * d_Up0(mdl.doffree);
                else
                    % Generalized Stiffness Parameter
                    GSP = n2 / (d_Up0(mdl.doffree)' * d_Up0_old(mdl.doffree));

                    % Adjust increment sign
                    if (GSP < 0)
                        sign = -sign;
                    end

                    % Adjustment factor of increment size
                    if (this.adjustStep == true)
                        J = sqrt(this.trg_iter/iter);
                    else
                        J = 1;
                    end

                    % Predicted increment of load ratio
                    d_lbd0 = this.method.predictedIncrement(this,mdl,sign,J,GSP,D_lbd,d_lbd0,D_U,d_Up0,Fref);
                end

                % Check increment of load ratio
                d_lbd0 = min(d_lbd0,this.max_increment);

                % Limit increment of load ratio to make total load ratio smaller than maximum value
                if ((this.max_lratio > 0.0 && lbd + d_lbd0 > this.max_lratio) ||...
                    (this.max_lratio < 0.0 && lbd + d_lbd0 < this.max_lratio))
                    d_lbd0 = this.max_lratio - lbd;
                end

                % Increments of load ratio and displacements for predicted solution
                d_lbd = d_lbd0;
                d_U0  = d_lbd0 * d_Up0;
                d_U   = d_U0;

                % Initialize incremental values of load ratio and displacements for current step
                D_lbd = d_lbd;
                D_U   = d_U;

                % Update total values of load ratio and displacements
                lbd = lbd + d_lbd;
                U   = U   + d_U;

                % Start iterative process
                iter = 1;
                conv = 0;

                while (conv == 0 && iter <= this.max_iter)
                    % Vector of external and internal forces
                    Fext = lbd * Fref;
                    [~,~,Fint,~,K] = mdl.globalMatrices(U);

                    % Vector of unbalanced forces
                    R = Fext - Fint;

                    % Check convergence
                    unbNorm = norm(R(mdl.doffree));
                    forNorm = norm(Fref(mdl.doffree));
                    conv = (unbNorm == 0 || forNorm == 0 || unbNorm/forNorm < this.tol);
                    if this.echo
                        fprintf(" iter.: %3d , ||R||/||F|| = %7.3e \n",iter,unbNorm/forNorm);
                    end
                    if conv == 1
                        break;
                    end

                    % Start/keep corrector phase
                    iter = iter + 1;

                    % Tangent and residual increments of displacements
                    d_Up = this.solveSystem(mdl,K,Fref);
                    d_Ur = this.solveSystem(mdl,K,R);

                    % Corrected increment of load ratio
                    d_lbd = this.method.correctedIncrement(this,mdl,d_lbd0,D_lbd,d_Up0_old,d_U0,d_Up,d_Ur,D_U,Fref,R);
                    if (~isreal(d_lbd))
                        conv = -1;
                        break;
                    end

                    % Corrected increment of displacements
                    d_U = d_lbd * d_Up + d_Ur;

                    % Increments of load ratio and displacements for current step
                    D_lbd = D_lbd + d_lbd;
                    D_U   = D_U   + d_U;

                    % Total values of load ratio and displacements
                    lbd = lbd + d_lbd;
                    U   = U   + d_U;
                end

                % Check for convergence fail or complex value of increment
                if (conv == 0)
                    if this.echo
                        disp('Convergence not achieved!');
                    end
                    return;
                elseif (conv == -1)
                    if this.echo
                        disp('Unable to compute load increment!');
                    end
                    return;
                end
                if this.echo
                    fprintf(' Step %d converged in iteration %-3d\n',step,iter);
                    fprintf(' Load factor: %f\n',lbd);
                end

                % Update state variables
                mdl.updateStateVar();

                % Store step results
                this.lbdplot(step+1) = lbd;
                this.Uplot(step+1) = U(mdl.ID(this.plotNd,this.plotDof));

                % Store predicted tangent increment of displacements for next step
                if (step ~= 1)
                    d_Up0_old = d_Up0;
                end

                % Check if maximum load ratio was reached
                if ((this.max_lratio >= 0 && lbd >= 0.999*this.max_lratio) ||...
                    (this.max_lratio <= 0 && lbd <= 0.999*this.max_lratio))
                    break;
                end
            end

            % Clean unused steps
            if (step < this.max_step)
                this.lbdplot = this.lbdplot(1:step+1);
                this.Uplot = this.Uplot(1:step+1);
            end

            % Save final result
            mdl.U = U;

            disp("*** Analysis completed!");
        end

        %------------------------------------------------------------------
        % Partition and solve a linear system of equations for free DOFs:
        %  f --> free DOF (natural B.C. - unknown) 
        %  c --> constrained DOF (essential B.C. - known) 
        %
        % [ Kff Kfs ] * [ Uf ] = [ Fext ]
        % [ Ksf Kss ]   [ Us ] = [   R  ]
        %
        function [U,Fext] = solveSystem(~,mdl,K,Fext,U)
            if nargin < 5
                U = zeros(mdl.ndof,1);
            end

            % Partition system of equations
            Kff = K(mdl.doffree,mdl.doffree);
            Ff  = Fext(mdl.doffree);

            % Solve system of equilibrium equations
            Uf = Kff \ Ff;

            % Displacement vector
            U(mdl.doffree) = Uf;
        end

        %------------------------------------------------------------------
        % Set the node and degree of freedom to be plotted.
        function setPlotDof(this,nd,dof)
            this.plotNd  = nd;
            this.plotDof = dof;
        end

        %------------------------------------------------------------------
        % Plot the load-displacement curve for the analysis.
        function plotCurves(this)
            figure;
            hold on, box on, grid on, axis on;
            plot(this.Uplot, this.lbdplot, 'o-k');
            xlabel('Displacement (m)');
            ylabel('Load factor');
            set(gca,'FontSize',16);
            set(gca,'FontName','Times');
        end
    end
end
