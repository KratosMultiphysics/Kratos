%% Model_HM Class
% This class defines a hydromechanical finite element model with 
% single-phase fluid flow. It extends the _Model_M_ class and incorporates 
% specific functionalities for hydromechanical analysis. Each node in the 
% model has three degrees of freedom:
% 
% * 2 displacement components (ux, uy)
% * 1 pore pressure (p)
% 
%% Methods
% * *setMaterial*: Sets the material properties using a _PorousMedia_ 
%                  object.
% * *initializeElements*: Initializes the elements of the model with their 
%                         properties..
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
% * *initializeDiscontinuitySegArray*: Initializes an array of 
%                                      discontinuity segments.
% * *initializeDiscontinuitySegment*: Initializes a single discontinuity 
%                                     segment.
% * *printResultsHeader*: Prints the header for the results table, 
%                         displaying node ID, displacement components 
%                         (ux, uy), and pore pressure (Pl).
% 
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class definition
classdef Model_HM < Model_M     
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        %% Embedded related data
        updateAperture = false;
        discontinuityTransversalFlow = false;
    end
    %% Constructor method
    methods
        function this = Model_HM()
            this = this@Model_M(false);
            this.ndof_nd = 3;       % Number of dofs per node
            this.physics = 'HM';    % Tag with the physics name
            this.condenseEnrDofs = false; % Enrichment dofs are global
            disp("*** Physics: Single-phase flow hydro-mechanical (HM)");
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Sets the material properties
        function setMaterial(this,porousMedia,fluid)
            if nargin < 3
                disp('Error in setMaterial: insufficient number of inputs.');
                disp('Physics HM requires 2 attribute(s): porousMedia, fluid.');
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
                udofs = this.getElementDofs(el,[1,2]);
                pdofs = this.getElementDofs(el,3);
                if (this.enriched == false)
                    elements(el) = RegularElement_HM(...
                                this.NODE(this.ELEM{el},:), this.ELEM{el},...
                                this.t, emat, this.intOrder,udofs,pdofs, ...
                                this.massLumping, this.lumpStrategy, this.isAxisSymmetric, ...
                                this.isPlaneStress);
                else
                    if (this.discontinuityTransversalFlow == true)
                        elements(el) = EnrichedElement_HM(...
                                this.NODE(this.ELEM{el},:), this.ELEM{el},...
                                this.t, emat, this.intOrder,udofs,pdofs, ...
                                this.massLumping, this.lumpStrategy, this.isAxisSymmetric, ...
                                this.isPlaneStress,this.addRelRotationMode, ...
                                this.addTangentialStretchingMode, this.addNormalStretchingMode,...
                                this.subDivIntegration, this.symmetricSDAEFEM);
                    else
                        elements(el) = EnrichedElementConductive_HM(...
                                this.NODE(this.ELEM{el},:), this.ELEM{el},...
                                this.t, emat, this.intOrder,udofs,pdofs, ...
                                this.massLumping, this.lumpStrategy, this.isAxisSymmetric, ...
                                this.isPlaneStress,this.addRelRotationMode, ...
                                this.addTangentialStretchingMode, this.addNormalStretchingMode,...
                                this.subDivIntegration, this.symmetricSDAEFEM);
                    end
                end
                if this.gravityOn
                    elements(el).type.gravityOn = true;
                end
            end
            this.element = elements;
        end

        % -----------------------------------------------------------------
        % Set the flag to update the aperture of the discontinuities
        function setUpdateAperture(this, flag)
            nDiscontinuities = this.getNumberOfDiscontinuities();
            for i = 1:nDiscontinuities
                nDiscontinuitySeg = this.discontinuitySet(i).getNumberOfDiscontinuitySegments();
                for j = 1:nDiscontinuitySeg
                    this.discontinuitySet(i).segment(j).updateAperture = flag;
                end
            end
        end

        % -----------------------------------------------------------------
        % Prescribe a pressure Dirichlet boundary condition at a node
        function resetPressureDirichletBC(this)
            this.resetDirichletBC(3);
        end

        % -----------------------------------------------------------------
        % Prescribe a pressure Dirichlet boundary condition at a node
        function setPressureDirichletBCAtDomain(this, value)
            this.setDirichletBCAtDomain(3, value);
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
        % Initializes an array of discontinuity segments
        function seg = initializeDiscontinuitySegArray(this,n)
            seg = [];
            if (this.discontinuityTransversalFlow == true)
                seg(n,1) = DiscontinuityElement_HM([],[]);
            else
                seg(n,1) = DiscontinuityElementConductive_HM([],[]);
            end
        end

        % -----------------------------------------------------------------
        % Initializes a single discontinuity segment
        function seg = initializeDiscontinuitySegment(this,nodeD,matD)
            if (this.discontinuityTransversalFlow == true)
                seg = DiscontinuityElement_HM(nodeD,matD);
            else
                seg = DiscontinuityElementConductive_HM(nodeD,matD);
            end
        end

    end
        %% Static methods
    methods (Static)

        % -----------------------------------------------------------------
        % Print header of the results
        function printResultsHeader()
            fprintf('\n  Node           ux        uy        Pl\n');
        end

    end
end
