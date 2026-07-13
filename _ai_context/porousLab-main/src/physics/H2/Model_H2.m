%% Model_H2 Class
% This class represents a two-phase flow finite element model.
% Each node has two degrees of freedom:
% 
% * 1 liquid phase pressure (Pl)
% * 1 gas phase pressure (Pg)
%
%% Methods
% * *setMaterial*: Sets the material properties using a _PorousMedia_ 
%                  object.
% * *initializeElements*: Initializes the elements of the model with their 
%                         properties.
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
%                         node ID, and pressures (Pl, Pg).
% 
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class definition
classdef Model_H2 < Model_H
    %% Constructor method
    methods
        function this = Model_H2()
            this = this@Model_H(false);
            this.ndof_nd = 2;       % Number of dofs per node
            this.physics = 'H2';    % Tag with the physics name
            disp("*** Physics: Two-phase hydraulic (H2)");
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Sets the material properties
        function setMaterial(this,porousMedia,liquidFluid,gasFluid)
            if nargin < 4
                disp('Error in setMaterial: insufficient number of inputs.');
                disp('Physics H2 requires 3 attribute(s): porousMedia, liquidFluid, gasFluid.');
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
            this.mat = struct('porousMedia',porousMedia,'liquidFluid',liquidFluid,'gasFluid',gasFluid);
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
                pl_dofs = this.getElementDofs(el,1);
                pg_dofs = this.getElementDofs(el,2);
                if (this.enriched == false)
                    elements(el) = RegularElement_H2(...
                            this.NODE(this.ELEM{el},:), this.ELEM{el},...
                            this.t, emat, this.intOrder, pl_dofs, pg_dofs, ...
                            this.massLumping, this.lumpStrategy, this.isAxisSymmetric);
                else
                    elements(el) = EnrichedElement_H2(...
                            this.NODE(this.ELEM{el},:), this.ELEM{el},...
                            this.t, emat, this.intOrder, pl_dofs, pg_dofs, ...
                            this.massLumping, this.lumpStrategy, this.isAxisSymmetric);
                end
                if this.gravityOn
                    elements(el).type.gravityOn = true;
                end
            end
            this.element = elements;
        end

        % -----------------------------------------------------------------
        % Prescribe a gas-phase pressure Dirichlet boundary condition at 
        % a node
        function setGasPressureDirichletBCAtNode(this, nodeId, value)
            this.setDirichletBCAtNode(nodeId, 2, value);
        end

        % -----------------------------------------------------------------
        % Prescribe a gas-phase pressure Dirichlet boundary condition at 
        % a point
        function setGasPressureDirichletBCAtPoint(this, X, value)
            this.setDirichletBCAtPoint(X, 2, value);
        end

        % -----------------------------------------------------------------
        % Prescribe a gas-phase pressure Dirichlet boundary condition at 
        % a border
        function setGasPressureDirichletBCAtBorder(this, border, value, range)
            if (nargin < 4), range = []; end
            this.setDirichletBCAtBorder(border, 2, value, range);
        end

        % -----------------------------------------------------------------
        % Prescribe a gas-phase pressure Neumann boundary condition at a 
        % node
        function setGasPressureNeumannBCAtNode(this, nodeId, value)
            this.setNeumannBCAtNode(nodeId, 2, value);
        end

        % -----------------------------------------------------------------
        % Prescribe a gas-phase pressure Neumann boundary condition at a 
        % point
        function setGasPressureNeumannBCAtPoint(this, X, value)
            this.setNeumannBCAtPoint(X, 2, value);
        end

        % -----------------------------------------------------------------
        % Prescribe a gas-phase pressure Neumann boundary condition at a 
        % border
        function setGasPressureNeumannBCAtBorder(this, border, value)
            this.setNeumannBCAtBorder(border, 2, value);
        end

        % -----------------------------------------------------------------
        % Sets the initial gas-phase pressure value for the whole domain
        function setInitialGasPressureAtDomain(this, value)
            this.setInitialDofAtDomain(2, value);
        end

        % -----------------------------------------------------------------
        % Sets the initial gas-phase pressure value at a node
        function setInitialGasPressureAtNode(this, nodeId, value)
            this.setInitialDofAtNode(nodeId, 2, value);
        end

        % -----------------------------------------------------------------
        % Initializes an array of discontinuity segments
        function seg = initializeDiscontinuitySegArray(~,n)
            seg(n,1) = DiscontinuityElement_H2([],[]);
        end

        % -----------------------------------------------------------------
        % Initializes a single discontinuity segment
        function seg = initializeDiscontinuitySegment(~,nodeD,matD)
            seg = DiscontinuityElement_H2(nodeD,matD);
        end
    end

    %% Static methods
    methods (Static)

        % -----------------------------------------------------------------
        % Print header of the results
        function printResultsHeader()
            fprintf('\n  Node           Pl        Pg\n');
        end

    end
end
