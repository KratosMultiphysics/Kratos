%% Model_M Class
% This class represents a mechanical finite element model. It extends the
% _Model_ class and is specifically designed for mechanical physics 
% simulations. Each node in the model has 2 degrees of freedom:
% 
% * 2 displacement components (ux, uy)
%
%% Methods
% * *setMaterial*: Sets the material properties using a _PorousMedia_ 
%                  object.
% * *initializeElements*: Initializes the elements of the model with their 
%                         properties.
% * *setDirichletBCAtNode*: Sets displacement Dirichlet boundary 
%                           conditions at a specific node.
% * *setDirichletBCAtPoint*: Sets displacement Dirichlet boundary 
%                            conditions at a specific point.
% * *setDirichletBCAtBorder*: Sets displacement Dirichlet boundary 
%                             conditions along a border.
% * *addLoadAtNode*: Adds a load at a specific node.
% * *addLoadAtPoint*: Adds a load at a specific point.
% * *addLoadAtBorder*: Adds a load along a specified border in a given 
%                      direction.
% * *initializeDiscontinuitySegArray*: Initializes an array of 
%                                      discontinuity segments.
% * *initializeDiscontinuitySegment*: Initializes a single discontinuity 
%                                     segment.
% * *addDiscontinuityData*: Adds additional discontinuity data.
% * *initializeDiscontinuitySegments*: Initializes discontinuity segments.
% * *plotDisplacementAlongSegment*: Plots displacement along a segment.
% * *plotDeformedMesh*: Plots the deformed mesh with a specified 
%                       amplification factor.
% * *updateResultVertices*: Updates the result vertices of each element 
%                           based on the specified configuration and 
%                           factor.
% * *printResultsHeader*: Prints the header for the results table.
% 
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class definition
classdef Model_M < Model   
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        %% Model data
        isPlaneStress       = false;
        %% Geostatic 
        addGeostatic        = false;
        addPorePressure     = false;
        stressfnc           = [];
        FeG                 = [];
        P                   = [];
        FeP                 = [];         
        %% Embedded related data
        addTangentialStretchingMode = false;
        addNormalStretchingMode     = false;
        addRelRotationMode          = false;
        symmetricSDAEFEM            = true;        % Flag for the use of a symmetric embedded formulation     
    end
    
    %% Constructor method
    methods
        function this = Model_M(printFlag)
            if nargin == 0, printFlag = true; end
            this = this@Model();
            this.ndof_nd = 2;       % Number of dofs per node
            this.physics = 'M';     % Tag with the physics name
            this.nNodalEnrDofs = 2;
            if (printFlag)
                disp("*** Physics: Mechanical");
            end 
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Sets the material properties
        function setMaterial(this,porousMedia)
            if nargin < 2
                disp('Error in setMaterial: insufficient number of inputs.');
                disp('Physics M requires 1 attribute(s): porousMedia.');
                error('Error in setMaterial.');
            end
            if ~isa(porousMedia,'PorousMedia')
                disp('Error in setMaterial: porousMedia is not a PorousMedia object.');
                error('Error in setMaterial.');
            end
            this.mat = struct('porousMedia',porousMedia);
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
                        'lc',this.getElementCharacteristicLength(el));
                dof_e = this.getElementDofs(el,[1,2]);
                if (this.enriched == false)
                    elements(el) = RegularElement_M(...
                                this.NODE(this.ELEM{el},:), this.ELEM{el},...
                                this.t, emat, this.intOrder,dof_e, ...
                                this.massLumping, this.lumpStrategy, this.isAxisSymmetric, ...
                                this.isPlaneStress);
                else
                    elements(el) = EnrichedElement_M(...
                                this.NODE(this.ELEM{el},:), this.ELEM{el},...
                                this.t, emat, this.intOrder,dof_e, ...
                                this.massLumping, this.lumpStrategy, this.isAxisSymmetric, ...
                                this.isPlaneStress,this.addRelRotationMode, ...
                                this.addTangentialStretchingMode, this.addNormalStretchingMode,...
                                this.condenseEnrDofs,this.subDivIntegration, this.symmetricSDAEFEM, ...
                                this.useNodalEnrDofs);
                end
                if this.gravityOn
                    elements(el).type.gravityOn = true;
                end
            end
            this.element = elements;
            
        end

        %------------------------------------------------------------------
        % Initialize additional model data associated with the physics
        function initializePhysicsAdditionalData(this)
            if this.addGeostatic
                this.initializeGeostaticStresses();
                this.initializeGeostaticForceVector();
            end
        end

        %------------------------------------------------------------------
        % Initialize the stresses in the elements with the given stress
        % state function
        function initializeGeostaticStresses(this)
            for el = 1 : this.nelem
                this.element(el).type.initialStress(this.stressfnc);
            end
        end

        %------------------------------------------------------------------
        % Initialize the vector with the geostatic forces
        function initializeGeostaticForceVector(this)
            [K,~,Fi,~] = this.globalMatrices(this.U);
            %this.FeG = sparse(this.ndof,1);
            this.FeG = Fi + K * this.U;
        end

        %------------------------------------------------------------------
        % Add the geostatic forces to the external force vector
        function Fe = addGeostaticForces(this, Fe)
            if ((this.addGeostatic == false) || isempty(this.FeG))
                return
            end
            Fe = Fe + this.FeG;
        end

        %------------------------------------------------------------------
        % Set the nodal pore-pressure values
        function setPorePressureField(this, P)
            this.P = P;
            if (isempty(this.FeP) == false)
                this.computePorePressureForceVct();
            end
        end

        %------------------------------------------------------------------
        % Add forces due to external pore-pressure field
        function Fe = addPorePressureForces(this, Fe)
            if ((this.addPorePressure == false) || (isempty(this.P)))
                return
            end
            if (isempty(this.FeP))
                this.computePorePressureForceVct();
            end
            Fe = Fe + this.FeP;
        end

        %------------------------------------------------------------------
        % Compute the forces due to the pore-pressure field.
        function computePorePressureForceVct(this)
            % Initialize
            this.FeP = zeros(this.ndof,1);
            % Compute and assemble element data
            for el = 1:this.nelem
                % Pressure "dofs"
                pdof = this.element(el).type.connect;
                % Displacement dofs
                udof = this.element(el).type.gle;
                % Get contribution from element el
                fep = this.element(el).type.porePressureForce(this.P(pdof));
                this.FeP(udof) = this.FeP(udof) + fep;
            end
            % Convert to sparse
            this.FeP = sparse(this.FeP);
        end

        %------------------------------------------------------------------
        % Add contribution of nodal loads to reference load vector.
        function Fe = addNodalLoad(this,Fe)
            Fe = addNodalLoad@Model(this,Fe);
            Fe = this.addGeostaticForces(Fe);
            Fe = this.addPorePressureForces(Fe);
        end

        % -----------------------------------------------------------------
        % Reset the displacements and strains of the model
        function resetDisplacements(this)
            for el = 1 : this.nelem
                udofs = this.element(el).type.resetDisplacements();
                this.U(udofs) = 0.0;
            end
        end

        % -----------------------------------------------------------------
        % Prescribe a displacement Dirichlet boundary condition at a node
        function setDisplacementDirichletBCAtNode(this, nodeId, value)
            this.setDirichletBCAtNode(nodeId, [1,2], value);
        end

        % -----------------------------------------------------------------
        % Prescribe a displacement Dirichlet boundary condition at a point
        function setDisplacementDirichletBCAtPoint(this, X, value)
            this.setDirichletBCAtPoint(X, [1,2], value);
        end

        % -----------------------------------------------------------------
        % Prescribe a displacement Dirichlet boundary condition at a border
        function setDisplacementDirichletBCAtBorder(this, border, value)
            this.setDirichletBCAtBorder(border, [1,2], value);
        end

        % -----------------------------------------------------------------
        % Impose a load at a node
        function addLoadAtNode(this, nodeId, value)
            this.setNeumannBCAtNode(nodeId, [1,2], value);
        end

        % -----------------------------------------------------------------
        % Impose a load at a point
        function addLoadAtPoint(this, X, value)
            this.setNeumannBCAtPoint(X, [1,2], value);
        end

        % -----------------------------------------------------------------
        % Impose a load at a border
        function addLoadAtBorder(this, border, dir, p, range)
            if nargin < 5
                if (strcmp(border,'left') || strcmp(border,'right')) 
                    range = [min(this.NODE(:,2)) , max(this.NODE(:,2))];
                elseif (strcmp(border,'top')||strcmp(border,'bottom'))
                    range = [min(this.NODE(:,1)) , max(this.NODE(:,1))];
                end
            end

            % Tolerance
            tol = 1.0e-12;

            % Get the nodes at the given border
            if strcmp(border,'left')
                ref = min(this.NODE(:,1));
                ndir = 1;
            elseif strcmp(border,'right')
                ref = max(this.NODE(:,1));
                ndir = 1;
            elseif strcmp(border,'top')
                ref = max(this.NODE(:,2));
                ndir = 2;
            elseif strcmp(border,'bottom')
                ref = min(this.NODE(:,2));
                ndir = 2;
            else
                disp('Warning: non-supported border.');
                disp('Available borders tag: ''left'',''right'', ''top'',''bottom''');
                return;
            end
            
            for el = 1:this.nelem 

                % Get number of linear interpolation points
                nLinNodes = length(this.ELEM{el});
                quadMesh  = false;
                if (nLinNodes == 6) || (nLinNodes == 8)
                    nLinNodes = nLinNodes / 2;
                    quadMesh  = true;
                end

                % Get the number of edges of the element
                nEdges = nLinNodes;
            
                % Get the coordinates of the element
                cX = [this.NODE(this.ELEM{el}(1:nLinNodes),1); this.NODE(this.ELEM{el}(1),1)];
                cY = [this.NODE(this.ELEM{el}(1:nLinNodes),2); this.NODE(this.ELEM{el}(1),2)];
            
                % Get the nodes of the borders
                NdBorders = [this.ELEM{el}(1:nLinNodes), this.ELEM{el}(1)];
            
                % Loop through the edges of the element
                for j = 1:nEdges
            
                    % coordinates of the edge
                    edgeX = [cX(j) , cX(j+1)];
                    edgeY = [cY(j) , cY(j+1)];
            
                    % select the edge
                    if ndir == 1
                        edge   = edgeX;
                        eRange = edgeY;
                    elseif ndir == 2
                        edge = edgeY;
                        eRange = edgeX;
                    end
            
                    % check if the edge belong to the boundary
                    if ((norm(edge-ref) < tol) && (min(eRange)>(range(1)-tol)) && (max(eRange)<(range(2)+tol)))
                        
                        % Compute the length of the edge
                        dx = edgeX(2) - edgeX(1);
                        dy = edgeY(2) - edgeY(1);
                        l = sqrt(dx*dx + dy*dy);
            
                        % id of the nodes of the edge
                        idNds = [NdBorders(j); NdBorders(j+1)];
            
                        % Equivalent nodal load
                        if quadMesh == false
                            feq = [0.5*p*l;0.5*p*l];
                        else
                            feq = [p*l;p*l;4.0*p*l]/6.0;
                            idNds = [idNds; this.ELEM{el}(j+nLinNodes)];
                        end

                        % For axisymmetric models
                        if this.isAxisSymmetric
                            feq = feq * 2 * pi * mean(edgeX);
                        end
            
                        % Add contribution to the LOAD matrix
                        this.LOAD(idNds,dir) = this.LOAD(idNds,dir) + feq;
                    end
                end
            end
        end

        % -----------------------------------------------------------------
        % Initializes an array of discontinuity segments
        function seg = initializeDiscontinuitySegArray(~,n)
            seg(n,1) = DiscontinuityElement_M([],[]);
        end

        % -----------------------------------------------------------------
        % Initializes a single discontinuity segment
        function seg = initializeDiscontinuitySegment(~,nodeD,matD)
            seg = DiscontinuityElement_M(nodeD,matD);
        end

        % -----------------------------------------------------------------
        % Adds some additional discontinuity data
        function addDiscontinuityData(this, additionalData)
            if nargin < 2 || ~isstruct(additionalData)
                return
            end
            fields = {
                'addRelRotationMode'
                'addTangentialStretchingMode'
                'addNormalStretchingMode'
            };
            for i = 1:numel(fields)
                f = fields{i};
                if isfield(additionalData, f) && ~isempty(additionalData.(f))
                    this.(f) = additionalData.(f);
                end
            end
        end

        %------------------------------------------------------------------
        % Initializes discontinuity segments
        function initializeDiscontinuitySegments(this)
            
            % First set the specific properties
            nDiscontinuities = this.getNumberOfDiscontinuities();
            for i = 1:nDiscontinuities
                nDiscontinuitySeg = this.discontinuitySet(i).getNumberOfDiscontinuitySegments();
                for j = 1:nDiscontinuitySeg
                    this.discontinuitySet(i).segment(j).addStretchingMode(this.addTangentialStretchingMode, this.addNormalStretchingMode);
                    this.discontinuitySet(i).segment(j).addRelRotationMode(this.addRelRotationMode);
                    this.discontinuitySet(i).segment(j).useNodalEnrDofs = this.useNodalEnrDofs;
                    if this.addPorePressure
                        this.discontinuitySet(i).segment(j).DP = 0.0;
                    end
                end
            end
            % Call the common initialization
            initializeDiscontinuitySegments@Model(this);
        end

        % -----------------------------------------------------------------
        % Plot the deformed mesh
        function plotDeformedMesh(this,amplFactor)
            this.updateResultVertices('Deformed',amplFactor);
            FEMPlot(this).plotMesh();
        end

        %------------------------------------------------------------------
        % Update the result nodes coordinates of each element
        function updateResultVertices(this,configuration,factor)
            for el = 1:this.nelem
                
                % Initialize the vertices array
                vertices = this.element(el).type.result.vertices0;

                % Get the updated vertices:
                if strcmp(configuration,'Deformed')

                    % Update the nodal displacement vector associated to the
                    % element. This displacement can contain the enhancement
                    % degrees of freedom.
                    this.element(el).type.ue = this.U(this.element(el).type.gle); 

                    % Update the vertices based on the displacement vector
                    % associated to the element
                    for i = 1:length(this.element(el).type.result.faces)
                        X = vertices(i,:);
                        u = this.element(el).type.displacementField(X);
                        vertices(i,:) = X + factor*u';
                    end
                end
                this.element(el).type.result.setVertices(vertices);
            end
        end

        %------------------------------------------------------------------
        % Evaluate a field at in a point inside an element
        % Assumes that the point is in the element domain
        function fieldValue = evaluateField(this, field, el, X)
            fieldValue = evaluateField@Model(this, field, el, X);
            if strcmp(field,'PressureExt')
                pdof = this.element(el).type.connect;
                fieldValue = this.element(el).type.pressureField(X,this.P(pdof));
            end
        end
    end
    %% Static methods
    methods (Static)

        % -----------------------------------------------------------------
        % Prints the header for the results table
        function printResultsHeader()
            fprintf('\n  Node           ux        uy\n');
        end

    end
end
