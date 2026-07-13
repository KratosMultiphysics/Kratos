%% Model_H2M Class
% This class represents a hydromechanical finite element model with 
% two-phase fluid flow. It extends the _Model_M_ class and is designed to 
% handle problems involving coupled solid deformation and fluid flow
% (liquid and gas phases) in porous media. Each node in the model has four
% degrees of freedom:
%
% * 2 displacement components (ux,uy)
% * 1 liquid-phase pore pressure (Pl)
% * 1 gas-phase pore pressure (Pg)
%
%% Methods
% * *setMaterial*: Sets the material properties using a _PorousMedia_ 
%                  object.
% * *initializeElements*: Initializes the elements of the model with their 
%                         properties.
% * *setPressureDirichletBCAtNode*: Sets pressure Dirichlet boundary 
%                                   conditions at a specific node.
% * *setPressureDirichletBCAtPoint*: Sets pressure Dirichlet boundary 
%                                    conditions at a specific point.
% * *setPressureDirichletBCAtBorder*: Sets pressure Dirichlet boundary 
%                                     conditions at a specific border.
% * *setPressureNeumannBCAtNode*: Sets a Neumann boundary condition for 
%                                 pressure at a specific node.
% * *setPressureNeumannBCAtPoint*: Sets a Neumann boundary condition for 
%                                  pressure at a specific point.
% * *setPressureNeumannBCAtBorder*: Sets a Neumann boundary condition for 
%                                   pressure along a border.
% * *setInitialPressureAtDomain*: Sets the initial pressure value for the 
%                                 entire domain.
% * *setInitialPressureAtNode*: Sets the initial pressure value at a 
%                               specific node.
% * *setGasPressureDirichletBCAtNode*: Sets pressure Dirichlet boundary 
%                                      conditions for gas-phase pore 
%                                      pressure at a specific node.
% * *setGasPressureDirichletBCAtPoint*: Sets pressure Dirichlet boundary 
%                                      conditions for gas-phase pore 
%                                      pressure at a specific point.
% * *setGasPressureDirichletBCAtBorder*: Sets pressure Dirichlet boundary 
%                                      conditions for gas-phase pore 
%                                      pressure at a specific border.
% * *setGasPressureNeumannBCAtNode*: Sets pressure Neumann boundary
%                                      conditions for gas-phase pore 
%                                      pressure at a specific node.
% * *setGasPressureNeumannBCAtPoint*: Sets pressure Neumann boundary
%                                      conditions for gas-phase pore 
%                                      pressure at a specific point.
% * *setGasPressureNeumannBCAtBorder*: Sets pressure Neumann boundary
%                                      conditions for gas-phase pore 
%                                      pressure at a specific border.
% * *setInitialGasPressureAtDomain*: Sets the initial gas-phase pressure 
%                                    value for the entire domain.
% * *setInitialGasPressureAtNode*: Sets the initial gas-phase pressure 
%                                  value at a specific node.
% * *printResultsHeader*: Prints a header for the results table, showing
%                         node ID, displacements (ux, uy), and pressures 
%                         (Pl, Pg).
% 
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class definition
classdef Model_H2M < Model_M    
    %% Constructor method
    methods
        function this = Model_H2M()
            this = this@Model_M(false);
            this.ndof_nd = 4;        % Number of dofs per node
            this.physics = 'H2M';    % Tag with the physics name
            disp("*** Physics: Two-phase flow hydro-mechanical (H2M)");
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Sets the material properties
        function setMaterial(this,porousMedia,liquidFluid,gasFluid)
            if nargin < 4
                disp('Error in setMaterial: insufficient number of inputs.');
                disp('Physics H2M requires 3 attribute(s): porousMedia, liquidFluid, gasFluid.');
                error('Error in setMaterial.');
            end
            if ~isa(porousMedia,'PorousMedia')
                disp('Error in setMaterial: porousMedia is not a PorousMedia object.');
                error('Error in setMaterial.');
            end
            if ~isa(liquidFluid,'Fluid')
                disp('Error in setMaterial: liquidFluid is not a Fluid object.');
                error('Error in setMaterial.');
            end
            if ~isa(gasFluid,'Fluid')
                disp('Error in setMaterial: gasFluid is not a Fluid object.');
                error('Error in setMaterial.');
            end
            this.mat  = struct('porousMedia',porousMedia,'liquidFluid',liquidFluid,'gasFluid',gasFluid);
        end

        %------------------------------------------------------------------
        % Initializes the elements of the model with the corresponding
        % properties
        function initializeElements(this)
            % Initialize the vector with the Element's objects
            elements(this.nelem,1) = Element(); 

            % Assemble the properties to the elements' objects
            for el = 1 : this.nelem
                % Create the material for the element
                emat =struct( ...
                        'porousMedia',this.mat.porousMedia(this.matID(el)), ...
                        'liquidFluid',this.mat.liquidFluid,...
                        'gasFluid',this.mat.gasFluid);
                u_dofs = this.getElementDofs(el,[1,2]);
                pl_dofs = this.getElementDofs(el,3);
                pg_dofs = this.getElementDofs(el,4);
                elements(el) = RegularElement_H2M(...
                            this.NODE(this.ELEM{el},:), this.ELEM{el},...
                            this.t, emat, this.intOrder,u_dofs,pl_dofs,pg_dofs, ...
                            this.massLumping, this.lumpStrategy, this.isAxisSymmetric, ...
                            this.isPlaneStress);
                if this.gravityOn
                    elements(el).type.gravityOn = true;
                end
            end
            this.element = elements;
        end

        % -----------------------------------------------------------------
        % Prescribe a pressure Dirichlet boundary condition at a node
        function setPressureDirichletBCAtNode(this, nodeId, value)
            this.setDirichletBCAtNode(nodeId, 3, value);
        end

        % -----------------------------------------------------------------
        % Prescribe a pressure Dirichlet boundary condition at a point
        function setPressureDirichletBCAtPoint(this, X, value)
            this.setDirichletBCAtPoint(X, 3, value);
        end

        % -----------------------------------------------------------------
        % Prescribe a pressure Dirichlet boundary condition at a border
        function setPressureDirichletBCAtBorder(this, border, value)
            this.setDirichletBCAtBorder(border, 3, value);
        end

        % -----------------------------------------------------------------
        % Prescribe a pressure Neumann boundary condition at a node
        function setPressureNeumannBCAtNode(this, nodeId, value)
            this.setNeumannBCAtNode(nodeId, 3, value);
        end

        % -----------------------------------------------------------------
        % Prescribe a pressure Neumann boundary condition at a point
        function setPressureNeumannBCAtPoint(this, X, value)
            this.setNeumannBCAtPoint(X, 3, value);
        end

        % -----------------------------------------------------------------
        % Prescribe a pressure Neumann boundary condition at a border
        function setPressureNeumannBCAtBorder(this, border, value)
            this.setNeumannBCAtBorder(border, 3, value);
        end

        % -----------------------------------------------------------------
        % Sets the initial pressure value for the whole domain
        function setInitialPressureAtDomain(this, value)
            this.setInitialDofAtDomain(3, value);
        end

        % -----------------------------------------------------------------
        % Sets the initial pressure value at a node
        function setInitialPressureAtNode(this, nodeId, value)
            this.setInitialDofAtNode(nodeId, 3, value);
        end

        % -----------------------------------------------------------------
        % Prescribe a gas-phase pressure Dirichlet boundary condition at 
        % a node
        function setGasPressureDirichletBCAtNode(this, nodeId, value)
            this.setDirichletBCAtNode(nodeId, 4, value);
        end

        % -----------------------------------------------------------------
        % Prescribe a gas-phase pressure Dirichlet boundary condition at 
        % a point
        function setGasPressureDirichletBCAtPoint(this, X, value)
            this.setDirichletBCAtPoint(X, 4, value);
        end

        % -----------------------------------------------------------------
        % Prescribe a gas-phase pressure Dirichlet boundary condition at 
        % a border
        function setGasPressureDirichletBCAtBorder(this, border, value)
            this.setDirichletBCAtBorder(border, 4, value);
        end

        % -----------------------------------------------------------------
        % Prescribe a gas-phase pressure Neumann boundary condition at a 
        % node
        function setGasPressureNeumannBCAtNode(this, nodeId, value)
            this.setNeumannBCAtNode(nodeId, 4, value);
        end

        % -----------------------------------------------------------------
        % Prescribe a gas-phase pressure Neumann boundary condition at a 
        % point
        function setGasPressureNeumannBCAtPoint(this, X, value)
            this.setNeumannBCAtPoint(X, 4, value);
        end

        % -----------------------------------------------------------------
        % Prescribe a gas-phase pressure Neumann boundary condition at a 
        % border
        function setGasPressureNeumannBCAtBorder(this, border, value)
            this.setNeumannBCAtBorder(border, 4, value);
        end

        % -----------------------------------------------------------------
        % Sets the initial gas-phase pressure value for the whole domain
        function setInitialGasPressureAtDomain(this, value)
            this.setInitialDofAtDomain(4, value);
        end

        % -----------------------------------------------------------------
        % Sets the initial gas-phase pressure value at a node
        function setInitialGasPressureAtNode(this, nodeId, value)
            this.setInitialDofAtNode(nodeId, 4, value);
        end
    end
        %% Static methods
    methods (Static)

        % -----------------------------------------------------------------
        % Print header of the results
        function printResultsHeader()
            fprintf('\n  Node           ux        uy        Pl        Pg\n');
        end

    end
end
