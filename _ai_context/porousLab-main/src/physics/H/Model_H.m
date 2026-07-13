%% Model_H Class
% This class represents a finite element model for single-phase hydraulic
% problems. Each node in the model has one degree of freedom:
%
% * 1 liquid phase pressure (Pl)
%
%% Methods
% * *setMaterial*: Sets the material properties using a _PorousMedia_ 
%                  object.
% * *initializeElements*: Initializes the elements of the model with their 
%                         properties.
% * *setPressureDirichletBCAtNode*: Sets pressure Dirichlet boundary 
%                                   conditions for liquid-phase pore 
%                                   pressure at a specific node.
% * *setPressureDirichletBCAtPoint*: Sets pressure Dirichlet boundary 
%                                    conditions for liquid-phase pore 
%                                    pressure at a specific point.
% * *setPressureDirichletBCAtBorder*: Sets pressure Dirichlet boundary 
%                                     conditions for liquid-phase pore 
%                                     pressure at a specific border.
% * *setPressureNeumannBCAtNode*: Sets pressure Neumann boundary
%                                 conditions for liquid-phase pore 
%                                 pressure at a specific node.
% * *setPressureNeumannBCAtPoint*: Sets pressure Neumann boundary
%                                  conditions for liquid-phase pore 
%                                  pressure at a specific point.
% * *setPressureNeumannBCAtBorder*: Sets pressure Neumann boundary
%                                   conditions for liquid-phase pore 
%                                   pressure at a specific border.
% * *setInitialPressureAtDomain*: Sets the initial pressure value for the 
%                                entire domain.
% * *setInitialPressureAtNode*: Sets the initial pressure value at a 
%                               specific node.
% * *initializeDiscontinuitySegArray*: Initializes an array of 
%                                      discontinuity elements with size n
%
% * *initializeDiscontinuitySegment*: Initializes a single discontinuity 
%                                     element with specified nodes and
%                                     material properties.
%
% * *createMaterialDataStructure*: Creates a data structure containing 
%                                  material properties such as fluid and 
%                                  initial aperture.
%
% * *printResultsHeader*: Prints a header for the results table, showing
%                         node ID and liquid pressure (Pl).
% 
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class definition
classdef Model_H < Model    
    properties (SetAccess = public, GetAccess = public)
        equivalentContinuum = false;
    end
    %% Constructor method
    methods
        function this = Model_H(printFlag)
            if nargin == 0, printFlag = true; end
            this = this@Model();
            this.ndof_nd = 1;             % Number of dofs per node
            this.physics = 'H';           % Tag with the physics name
            this.condenseEnrDofs = false; % Enrichment dofs are global
            if (printFlag)
                disp("*** Physics: Single-phase hydraulic (H)");
            end
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Sets the material properties
        function setMaterial(this,porousMedia,fluid)
            if nargin < 3
                disp('Error in setMaterial: insufficient number of inputs.');
                disp('Physics H requires 2 attribute(s): porousMedia, fluid.');
                error('Error in setMaterial.');
            end
            if ~isa(porousMedia,'PorousMedia')
                disp('Error in setMaterial: porousMedia is not a PorousMedia object.');
                error('Error in setMaterial.');
            end
            if ~isa(fluid,'Fluid')
                disp('Error in setMaterial: fluid is not a Fluid object.');
                error('Error in setMaterial.');
            end
            this.mat = struct('porousMedia',porousMedia,'fluid',fluid);
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
                        'fluid',this.mat.fluid);
                dof_e = this.getElementDofs(el,1);
                if (this.enriched == false)
                    elements(el) = RegularElement_H(...
                                this.NODE(this.ELEM{el},:), this.ELEM{el},...
                                this.t, emat, this.intOrder,dof_e, ...
                                this.massLumping, this.lumpStrategy, this.isAxisSymmetric);
                else
                    elements(el) = EnrichedElement_H(...
                                    this.NODE(this.ELEM{el},:), this.ELEM{el},...
                                    this.t, emat, this.intOrder,dof_e, ...
                                    this.massLumping, this.lumpStrategy, this.isAxisSymmetric);
                end
                if this.gravityOn
                    elements(el).type.gravityOn = true;
                end
            end
            this.element = elements;
        end

        %------------------------------------------------------------------
        % Initialize additional model data associated with the physics
        function initializePhysicsAdditionalData(~)
        end

        % -----------------------------------------------------------------
        % Prescribe a pressure Dirichlet boundary condition at a node
        function setPressureDirichletBCAtNode(this, nodeId, value)
            this.setDirichletBCAtNode(nodeId, 1, value);
        end

        % -----------------------------------------------------------------
        % Prescribe a pressure Dirichlet boundary condition at a point
        function setPressureDirichletBCAtPoint(this, X, value)
            this.setDirichletBCAtPoint(X, 1, value);
        end

        % -----------------------------------------------------------------
        % Prescribe a pressure Dirichlet boundary condition at a border
        function setPressureDirichletBCAtBorder(this, border, value, range)
            if (nargin < 4), range = []; end
            this.setDirichletBCAtBorder(border, 1, value, range);
        end

        % -----------------------------------------------------------------
        % Prescribe a pressure Neumann boundary condition at a node
        function setPressureNeumannBCAtNode(this, nodeId, value)
            this.setNeumannBCAtNode(nodeId, 1, value);
        end

        % -----------------------------------------------------------------
        % Prescribe a pressure Neumann boundary condition at a point
        function setPressureNeumannBCAtPoint(this, X, value)
            this.setNeumannBCAtPoint(X, 1, value);
        end

        % -----------------------------------------------------------------
        % Prescribe a pressure Neumann boundary condition at a border
        function setPressureNeumannBCAtBorder(this, border, value, range)
            if (nargin < 4), range = []; end
            this.setNeumannBCAtBorder(border, 1, value, range);
        end

        % -----------------------------------------------------------------
        % Sets the initial pressure value for the whole domain
        function setInitialPressureAtDomain(this, value)
            this.setInitialDofAtDomain(1, value);
        end

        % -----------------------------------------------------------------
        % Sets the initial pressure value at a node
        function setInitialPressureAtNode(this, nodeId, value)
            this.setInitialDofAtNode(nodeId, 1, value);
        end

        % -----------------------------------------------------------------
        % Initializes an array of discontinuity segments
        function seg = initializeDiscontinuitySegArray(~,n)
            seg(n,1) = DiscontinuityElement_H([],[]);
        end

        % -----------------------------------------------------------------
        % Initializes a single discontinuity segment
        function seg = initializeDiscontinuitySegment(~,nodeD,matD)
            seg = DiscontinuityElement_H(nodeD,matD);
        end

        %------------------------------------------------------------------
        % Create material data structure
        function mat = createMaterialDataStructure(this)
            mat = struct( ...
                'fluid',this.fluid, ...
                'initialAperture',this.initialAperture);
        end
    end

    %% Static methods
    methods (Static)

        % -----------------------------------------------------------------
        % Print header of the results
        function printResultsHeader()
            fprintf('\n  Node           Pl\n');
        end

    end
end
