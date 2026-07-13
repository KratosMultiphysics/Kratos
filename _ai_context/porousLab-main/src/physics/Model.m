%% Model Class
% Abstract class to create and manage a Finite Element model. It provides 
% methods and properties to define the mesh, boundary conditions, material 
% properties, and other essential components of a FEM simulation. 
% The class also includes methods for pre-processing, assembling global 
% matrices, applying boundary conditions, and updating state variables.
%
%% Methods
% * *setMaterial*: Sets the material object.
% * *initializeElements*: Initializes the element objects.
% * *printResultsHeader*: Configures the header for printing results.
% * *setMesh*: Sets the mesh nodes and connectivity.
% * *initializeBasicVariables*: Initializes basic variables like number of 
%                               nodes, elements, and DOFs.
% * *checkMaterialId*: Ensures the material ID vector is initialized.
% * *createNodeDofIdMtrx*: Creates the ID matrix and determines free and 
%                          fixed DOFs.
% * *setDirichletBCAtNode*: Sets Dirichlet boundary conditions at a 
%                           specific node.
% * *setDirichletBCAtPoint*: Sets Dirichlet boundary conditions at a 
%                            specific point.
% * *setDirichletBCAtBorder*: Sets Dirichlet boundary conditions along a 
%                             border.
% * *setNeumannBCAtNode*: Sets Neumann boundary conditions at a specific 
%                         node.
% * *setNeumannBCAtPoint*: Sets Neumann boundary conditions at a specific 
%                          point.
% * *setNeumannBCAtBorder*: Sets Neumann boundary conditions along a 
%                           border.
% * *setInitialDofAtDomain*: Sets initial DOF values for the entire domain.
% * *setInitialDofAtNode*: Sets initial DOF values at a specific node.
% * *getNodesAtBorder*: Retrieves nodes along a specified border.
% * *closestNodeToPoint*: Finds the closest node to a given point.
% * *preComputations*: Performs pre-processing tasks like initializing 
%                      elements and assembling matrices.
% * *getElementDofs*: Retrieves DOFs for a specific element.
% * *initializeDisplacementVct*: Initializes the global displacement 
%                                vector.
% * *assembleDiscontinuitySegments*: Assembles discontinuity segments into 
%                                    elements.
% * *addNodalLoad*: Adds nodal loads to the reference load vector.
% * *getElementsCharacteristicLength*: Obtain the characteristic length of
%                                      the finite element.
% * *getElementCharacteristicLength*: Computes the characteristic length 
%                                     for a specific element.
% * *getNodeCharacteristicLength*: Compute the mean characteristic length
%                                  of the elements associated with each
%                                  node.
% * *initializeSparseMtrxAssemblageVariables*: Initializes variables for 
%                                              sparse matrix assembly.
% * *globalMatrices*: Assembles global system matrices and vectors.
% * *getLinearSystem*: Assembles the linear system for solving.
% * *applyDirichletBC*: Applies Dirichlet boundary conditions to the 
%                       system.
% * *updateStateVar*: Updates state variables at integration points.
% * *resequenceNodes*: Resequences nodes using the reverse Cuthill-McKee 
%                      algorithm.
% * *rebuildConnectivity*: Rebuilds connectivity matrices after 
%                          resequencing.
% * *addPreExistingDiscontinuities*: Adds pre-existing discontinuities to 
%                                    the model.
% * *initializeDiscontinuitySegments*: Initializes discontinuity segments.
% * *getNumberOfDiscontinuities*: Returns the number of discontinuities.
% * *useEnrichedFormulation*: Enables or disables enriched formulation.
% * *printResults*: Prints nodal results.
% * *evaluateField*: Evaluate a field at in a point inside an element. It
%                    assumes that the point is in the element domain.
% * *updateResultVertexData*: Updates result data for each element's 
%                             vertices.
% * *plotField*: Plots a specified field over the mesh.
% * *plotFieldAlongSegment*: Plot a given field along a given segment.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00: Initial version (April 2023).
% 
%% Class definition
classdef Model < handle    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        name                = 'mdl';         % Model name
        physics             = [];            % Physics of the problem
        NODE                = [];            % Nodes of the fem mesh
        ELEM                = [];            % Nodes connectivity
        DIRICHLET_TAG       = [];            % Matrix with tags indicating the Dirichlet BC
        DIRICHLET_VAL       = [];            % Matrix with values of the Dirichlet BC
        LOAD                = [];            % Matrix with the nodal Neumann BC
        INIT                = [];            % Matrix with the initial values of each dof
        t                   = 1.0;           % Thickness
        mat                 = [];            % Struct with material properties
        intOrder            = 2;             % Order of the numerical integration quadrature
        nnodes              = 1;             % Number of nodes
        nelem               = 1;             % Number of elements
        doffree             = [];            % Vector with the free dofs
        doffixed            = [];            % Vector with the fixed dofs
        ndof_nd             = 2;             % Number of dof per node
        ndof                = 1;             % Number of degrees of freedom
        ndoffree            = 0;             % Number of free degrees of freedom
        ndoffixed           = 0;             % Number of fixed degrees of freedom
        Dof                 = [];            % Vector with all regular dofs
        ID                  = [];            % Each line of the ID matrix contains the global numbers for the node DOFs
        U                   = [];            % Global displacement vector
        element             = [];            % Array with the element's objects
        nDofElemTot         = 0;             % Aux value used to sparse matrix assemblage
        sqrNDofElemTot      = 0;             % Aux value used to sparse matrix assemblage
        matID               = [];            % Vector with the material id of each element
        gravityOn           = false;         % Flag to consider the gravity forces
        massLumping         = false;         % Tag for applying a mass lumping process
        lumpStrategy        = 1;             % Id of the mass lumping strategy
        isAxisSymmetric     = false;         % Flag for axissymetric models
        enriched            = false;         % Flag to use embedded formulation
        discontinuitySet    = [];            % Array with the discontinuity objects
        condenseEnrDofs     = true;          % Flag to condense the enrichment dofs
        dofenr              = [];            % Vector with the enrichment dofs
        ndofenr             = 0;             % Number of enrichment dofs
        useNodalEnrDofs     = false;         % Flag to use nodal enrichment dofs
        nNodalEnrDofs       = 0;             % Number of nodal enrichment dofs (used if useNodalEnrDofs == true)
        subDivIntegration   = false;         % Flag to apply a sub-division of the element to define the integration points
        initializeMdl       = false;         % Flag to check if the model has been initialized
    end
    
    %% Constructor method
    methods
        function this = Model()
            disp("*** Initializing model...")
        end
    end

    %% Abstract methods
    methods(Abstract)

        % Set the material object
        setMaterial(this, varargin)

        % Initialize the elements objects
        initializeElements(this);

        % Initialize additional model data associated with the physics
        initializePhysicsAdditionalData(this);

        % Configure the header to printed when printing results
        printResultsHeader();
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        function setMesh(this,node,elem)
            % Set the mesh nodes coordinates and connectivity
            this.NODE = node;
            this.ELEM = elem;

            % Initialize basic variables
            this.initializeBasicVariables();

        end

        %------------------------------------------------------------------
        % Initializes the basic variables of the model
        function initializeBasicVariables(this)
            this.nnodes        = size(this.NODE,1);
            this.nelem         = size(this.ELEM,1);           
            this.ndof          = this.ndof_nd * this.nnodes; 
            this.DIRICHLET_TAG = zeros(this.nnodes,this.ndof_nd);
            this.DIRICHLET_VAL = zeros(this.nnodes,this.ndof_nd);
            this.LOAD          = zeros(this.nnodes,this.ndof_nd);
            this.INIT          = zeros(this.nnodes,this.ndof_nd);
        end

        %------------------------------------------------------------------
        % Check is the material is well defined or not
        function checkMaterialId(this)
            if isempty(this.matID)
                this.matID  = ones(this.nelem,1);
            end
        end

        %------------------------------------------------------------------
        % Creates and assembles the matrix containing all the degrees of
        % freedom
        function createNodeDofIdMtrx(this)
            % Initialize the ID matrix and the number of fixed dof
            this.ID = zeros(this.nnodes,this.ndof_nd);
            this.ndoffixed = 0;
            
            % Assemble the ID matrix
            for i = 1:this.nnodes
                for j = 1:this.ndof_nd
                    this.ID(i,j) = (i - 1) * this.ndof_nd + j;
                    if (this.DIRICHLET_TAG(i,j) == 1)
                        this.ndoffixed = this.ndoffixed + 1;
                    end
                end
            end

            % Vector with all the dofs
            this.Dof = 1:this.ndof;
            
            % Number of free dof
            this.ndoffree = this.ndof - this.ndoffixed;
            
            % Initialize the counters
            this.doffixed = zeros(this.ndoffixed,1);
            this.doffree  = zeros(this.ndoffree,1);
            
            % Update the ID matrix with the free dof numbered first
            countFree = 1;
            countFixed = 1;
            for i = 1:this.nnodes
                for j = 1:this.ndof_nd
                    if this.DIRICHLET_TAG(i,j) == 1
                        this.doffixed(countFixed) = this.ID(i,j);
                        countFixed = countFixed + 1;
                    else 
                        this.doffree(countFree) = this.ID(i,j);
                        countFree = countFree + 1;
                    end
                end
            end

            % Add the enrichment dofs to the free dof vector
            for i = 1:this.ndofenr
                this.doffree(countFree) = this.dofenr(i);
                countFree = countFree + 1;
            end
        end

        %------------------------------------------------------------------
        % Prescribe a Dirichlet boundary condition at a node        
        function resetDirichletBC(this, dofId)
            for i = 1:this.nnodes
                for j = 1:length(dofId)
                    this.DIRICHLET_TAG(i,dofId(j)) = 0;
                    this.DIRICHLET_VAL(i,dofId(j)) = NaN;
                end
            end
        end
        
        %------------------------------------------------------------------
        % Prescribe a Dirichlet boundary condition at a node
        function setDirichletBCAtNode(this, nodeId, dofId, value)
            if (length(dofId) ~= length(value))
                disp('Error prescribing Dirichlet BC at a node');
                disp('length(dofId) ~= length(value)');
                error('Error in setDirichletBCAtNode');
            end
            for i = 1:length(dofId)
                if ( this.DIRICHLET_TAG(nodeId,dofId(i)) == 0)
                    this.DIRICHLET_TAG(nodeId,dofId(i)) = ~isnan(value(i));
                    this.DIRICHLET_VAL(nodeId,dofId(i)) = value(i);
                else
                    if ((this.DIRICHLET_VAL(nodeId,dofId(i)) ~= value(i)) && ~isnan(value(i)))
                        disp([' ** Warning: this node had a different prescribed value:',num2str(nodeId)])
                        this.DIRICHLET_VAL(nodeId,dofId(i)) = value(i);
                    end
                end
            end
        end

        %------------------------------------------------------------------
        % Prescribe a Dirichlet boundary condition at a node        
        function setDirichletBCAtDomain(this, dofId, value)
            for i = 1:this.nnodes
                this.setDirichletBCAtNode(i,dofId,value);
            end
        end

        %------------------------------------------------------------------
        % Prescribe a Dirichlet boundary condition at a node        
        function setDirichletBCAtPoint(this, X, dofId, value)
            nodeId = this.closestNodeToPoint(X);
            this.setDirichletBCAtNode(nodeId,dofId,value);
        end

        %------------------------------------------------------------------
        % Prescribe a Dirichlet boundary condition at a border
        function setDirichletBCAtBorder(this, border, dofId, value, range)
            if ((nargin < 5) || isempty(range))
                if strcmp(border,'left') || strcmp(border,'right')
                    range = [min(this.NODE(:,2)) , max(this.NODE(:,2))];
                elseif strcmp(border,'top') || strcmp(border,'bottom')
                    range = [min(this.NODE(:,1)) , max(this.NODE(:,1))];
                end
            end
            nodeIds = this.getNodesAtBorder(border,range);
            for i = 1:length(nodeIds)
                this.setDirichletBCAtNode(nodeIds(i),dofId,value);
            end
        end

        %------------------------------------------------------------------
        % Prescribe a Neumann boundary condition at a node
        function setNeumannBCAtNode(this, nodeId, dofId, value)
            this.LOAD(nodeId,dofId) = value;
        end

        %------------------------------------------------------------------
        % Prescribe a Dirichlet boundary condition at a point
        function setNeumannBCAtPoint(this, X, dofId, value)
            nodeId = this.closestNodeToPoint(X);
            this.LOAD(nodeId,dofId) = value;
        end

        %------------------------------------------------------------------
        % Prescribe a Neumann boundary condition at a node
        function setNeumannBCAtBorder(this, border, dofId, value, range)
            if ((nargin < 5) || isempty(range))
                if strcmp(border,'left') || strcmp(border,'right')
                    range = [min(this.NODE(:,2)) , max(this.NODE(:,2))];
                elseif strcmp(border,'top') || strcmp(border,'bottom')
                    range = [min(this.NODE(:,1)) , max(this.NODE(:,1))];
                end
            end
            nodeIds = this.getNodesAtBorder(border, range);
            for i = 1:length(nodeIds)
                this.setNeumannBCAtNode(nodeIds(i),dofId,value);
            end
        end

        %------------------------------------------------------------------
        % Prescribe an initial boundary condition at the whole
        % domain
        function setInitialDofAtDomain(this, dofId, value)
            if (length(dofId) ~= length(value))
                disp('Error setting initial dof value at the domain');
                disp('length(dofId) ~= length(value)');
                error('Error in setInitialDofAtDomain');
            end
            this.INIT(:,dofId) = value;
        end

        %------------------------------------------------------------------
        % Prescribe an initial boundary condition at a node
        function setInitialDofAtNode(this, nodeId, dofId, value)
            if (length(dofId) ~= length(value))
                disp('Error setting initial dof value at the domain');
                disp('length(dofId) ~= length(value)');
                error('Error in setInitialDofAtNode');
            end
            this.INIT(nodeId,dofId) = value;
        end

        %------------------------------------------------------------------
        % Identify the nodes contained in any of the borders
        function nodeIds = getNodesAtBorder(this,border,range)
            if ((nargin < 3) || isempty(range))
                if strcmp(border,'left') || strcmp(border,'right')
                    range = [min(this.NODE(:,2)) , max(this.NODE(:,2))];
                elseif strcmp(border,'top') || strcmp(border,'bottom')
                    range = [min(this.NODE(:,1)) , max(this.NODE(:,1))];
                end
            end
            % Get the nodes at the given border
            if strcmp(border,'left')
                nodeIds = find((abs(this.NODE(:,1)-min(this.NODE(:,1)))<1.0e-12) & ((this.NODE(:,2))>range(1)-1.0e-12) & ((this.NODE(:,2))<range(2)+1.0e-12));
            elseif strcmp(border,'right')
                nodeIds = find((abs(this.NODE(:,1)-max(this.NODE(:,1)))<1.0e-12) & ((this.NODE(:,2))>range(1)-1.0e-12) & ((this.NODE(:,2))<range(2)+1.0e-12));
            elseif strcmp(border,'top')
                nodeIds = find((abs(this.NODE(:,2)-max(this.NODE(:,2)))<1.0e-12) & ((this.NODE(:,1))>range(1)-1.0e-12) & ((this.NODE(:,1))<range(2)+1.0e-12));
            elseif strcmp(border,'bottom')
                nodeIds = find((abs(this.NODE(:,2)-min(this.NODE(:,2)))<1.0e-12) & ((this.NODE(:,1))>range(1)-1.0e-12) & ((this.NODE(:,1))<range(2)+1.0e-12));
            else
                disp('Warning: non-supported border.');
                disp('Available borders tag: ''left'',''right'', ''top'',''bottom''');
                nodeIds = [];
            end
        end

        %------------------------------------------------------------------
        % Update a prescribed Dirichlet boundary condition value at a node
        function updateValueDirichletBCAtNode(this, nodeId, dofId, value)
            if (length(dofId) ~= length(value))
                disp('Error updating prescribed Dirichlet BC at a node');
                disp('length(dofId) ~= length(value)');
                error('Error in updateDirichletBCAtNode');
            end
            for i = 1:length(dofId)
                if ( this.DIRICHLET_TAG(nodeId,dofId(i)) == 1)
                    this.DIRICHLET_VAL(nodeId,dofId(i)) = value(i);
                else
                    disp('Error updating prescribed Dirichlet BC at a node');
                    disp('This node did not have a prescribed value.')
                    error('Error in updateDirichletBCAtNode');
                end
            end
        end

        %------------------------------------------------------------------
        % Update the DOFs vector to respect the new BC values.
        function updateDirichletBC(this)
            this.applyDirichletBCtoDOFVct();
            this.createNodeDofIdMtrx();
        end
        
        %------------------------------------------------------------------
        % Find the closest node to a given point
        function nd = closestNodeToPoint(this,X)
            if size(X,1) == 2, X = X'; end
            d = vecnorm((this.NODE - X)');
            [~,id] = sort(d);
            nd = id(1);
        end

        %------------------------------------------------------------------
        % Perform all the necessary pre-computations to initialize the 
        % model
        function preComputations(this)
            if(this.initializeMdl == false)
                disp("*** Pre-processing...");
                
                % Check and initialize the material ID vector
                this.checkMaterialId();
    
                % Create nodes DOF ids matrix
                this.createNodeDofIdMtrx();
    
                % Initialize elements
                this.initializeElements();
    
                % Assemble discontinuity segments to the elements
                this.assembleDiscontinuitySegments();

                % Initialize integration points
                this.initializeIntegrationPoints();
                
                % Compute auxiliar variables for assemblage of sparse matrices
                this.initializeSparseMtrxAssemblageVariables();
    
                % Initialize the displacement vector
                this.initializeDisplacementVct();

                % Initialize physics additional data
                this.initializePhysicsAdditionalData();
    
                % Update flag to indicate that the model has already been initialized
                this.initializeMdl = true;
            end
        end

        %------------------------------------------------------------------
        % Initialize the elements integration points
        function initializeIntegrationPoints(this)
            for el = 1 : this.nelem
                this.element(el).type.initializeIntPoints();
            end
        end

        %------------------------------------------------------------------
        % Reset the model dof vector and state variables
        function resetModelState(this)
            if(this.initializeMdl == true)
                % Reset the displacement vector
                this.initializeDisplacementVct();
                % Reset the state variables
                for el = 1:this.nelem
                    this.element(el).type.resetIntegrationPts();
                end
            end
        end

        %------------------------------------------------------------------
        % Obtain all the element degrees of freedom
        function dof = getElementDofs(this,el,dofId)
            nnd_el = length(this.ELEM{el});
            dof = reshape(this.ID(this.ELEM{el},dofId)', 1, nnd_el*length(dofId));
        end

        %------------------------------------------------------------------
        % Initialize the vector containing all the displacement
        % values
        function initializeDisplacementVct(this)
            % Initialize the displacement vector 
            this.U = zeros(this.ndof,1);

            % Set the initial values
            this.applyICtoDOFVct();

            % Set the prescribed values
            this.applyDirichletBCtoDOFVct();

            % Save initial dofs to the elements
            for el = 1 : this.nelem
                this.element(el).type.ue = this.U(this.element(el).type.gle);
            end
        end

        %------------------------------------------------------------------
        % Apply the initial conditions to the DOFs vector
        function applyICtoDOFVct(this)
            for i = 1:this.nnodes
                for j = 1:this.ndof_nd
                    this.U(this.ID(i,j)) = this.INIT(i,j);
                end
            end
        end

        %------------------------------------------------------------------
        % Apply the Dirichlet BC to the DOFs vector
        function applyDirichletBCtoDOFVct(this)
            for i = 1:this.nnodes
                for j = 1:this.ndof_nd
                    if (this.DIRICHLET_TAG(i,j) == 1.0)
                        this.U(this.ID(i,j)) = this.DIRICHLET_VAL(i,j);
                    end
                end
            end
        end

        %------------------------------------------------------------------
        % Assembly discontinuity segments for enriched elements
        function assembleDiscontinuitySegments(this)
            if (this.enriched == false)
                return
            elseif isempty(this.discontinuitySet)
                return
            else
                nDiscontinuities = length(this.discontinuitySet);
                for i = 1:nDiscontinuities
                    % Loop through the segments of this discontinuity
                    k = 1;
                    for j = 1:size(this.discontinuitySet(i).Xlin,1)-1
                        el = this.discontinuitySet(i).elemID(j);
                        if (el > 0)
                            dseg = this.discontinuitySet(i).segment(k);
                            this.element(el).type.addDiscontinuitySegment(dseg);
                            k = k + 1;
                        end
                    end
                end
                this.updateEnrichedElementDofVector();
            end
        end
        
        %------------------------------------------------------------------
        % Add enrichment dofs to the elements dof vector
        function updateEnrichedElementDofVector(this)
            if (this.condenseEnrDofs)
                return
            else
                for el = 1:this.nelem
                    this.element(el).type.addEnrichmentToDofVector();
                end
            end
        end

        %------------------------------------------------------------------
        % Add contribution of nodal loads to reference load vector.
        function Fe = addNodalLoad(this,Fe)
            for i = 1:this.nnodes
                for j = 1:this.ndof_nd
                    Fe(this.ID(i,j)) = Fe(this.ID(i,j)) + this.LOAD(i,j);
                end
            end
        end

        %------------------------------------------------------------------
        % Obtain the characteristic length of the finite element
        function Lce = getElementsCharacteristicLength(this)
            Lce=zeros(this.nelem,1);
            for el = 1:this.nelem
                Lce(el) = this.getElementCharacteristicLength(el);
            end
        end

        %------------------------------------------------------------------
        % Get characteristic length of the elements
        function Lce = getElementCharacteristicLength(this,el)
            % Vertices of the element el coordinates
            vx = this.NODE(this.ELEM{el},1); 
            vy = this.NODE(this.ELEM{el},2);
        
            % Number of vertices 
            nv = length(this.ELEM{el}); 
        
            % Shifted vertices
            vxS = vx([2:nv 1]);
            vyS = vy([2:nv 1]); 
        
            % Compute the area of the element (trapezoidal rule)
            temp = vx.*vyS - vy.*vxS;
            Ae   = 0.5*sum(temp);
            
            % Characteristic lenght (quadrilateral elements)
            Lce = sqrt(Ae);
            if (nv == 3)||(nv == 6)
                Lce = Lce * sqrt(2.0);
            end
        end
        
        %------------------------------------------------------------------
        % Compute the mean characteristic length of the elements associated 
        % with each node
        function Lc = getNodeCharacteristicLength(this)
            Lce = this.getElementsCharacteristicLength();
            Lc = zeros(this.nnodes,1);
            for i = 1:this.nnodes
                % Get the elements associated with this node
                idElem = cellfun(@(elem) any(elem == i), this.ELEM);
                % Compute the mean characteristic lenght of these nodes
                Lc(i) = mean(Lce(idElem));
            end

        end

        %------------------------------------------------------------------
        % Initializes variables related to the assembly of sparse matrices 
        % in the model.
        function initializeSparseMtrxAssemblageVariables(this)
            this.nDofElemTot = 0;
            this.sqrNDofElemTot = 0;
            for el = 1:this.nelem
                this.nDofElemTot = this.nDofElemTot + this.element(el).type.ngle;
                this.sqrNDofElemTot = this.sqrNDofElemTot + this.element(el).type.ngle*this.element(el).type.ngle;
            end
        end

        %------------------------------------------------------------------
        % Global system matrices
        function [K, C, Fi, Fe, dfidx] = globalMatrices(this,U)   
            % Indices for the assembling the matrices and vector
            iDof = zeros(this.sqrNDofElemTot,1);
            jDof = zeros(this.sqrNDofElemTot,1);
            eDof = zeros(this.nDofElemTot,1);

            % Initialize the components of the matrices and vector
            K_ij      = zeros(this.sqrNDofElemTot,1);
            C_ij      = zeros(this.sqrNDofElemTot,1);
            dfidx_ij  = zeros(this.sqrNDofElemTot,1);
            fi_i      = zeros(this.nDofElemTot,1);
            fe_i      = zeros(this.nDofElemTot,1);

            % Initialize auxiliar variables
            counterK = 0;
            counterF = 0;
            
            % Compute and assemble element data
            for el = 1:this.nelem
                
                % Get the element object
                elementType = this.element(el).type;

                % Update the element displacement vector
                elementType.ue = U(elementType.gle);

                % Get the vector of the element dof  
                gle_i  = elementType.gle;
                gle_j  = gle_i;
                nGlei  = length(gle_i);
                nGlej  = length(gle_j);
                nGleij = nGlei*nGlej;

                % Get the indices for assemblage
                iDofEl = repmat(gle_i',1,nGlej);
                jDofEl = repmat(gle_j,nGlei,1);
                iDof(counterK+1:counterK+nGleij) = iDofEl(:);
                jDof(counterK+1:counterK+nGleij) = jDofEl(:);
                eDof(counterF+1:counterF+nGlei)  = gle_i';
            
                % Get local matrices and vector
                [K_e,C_e,fi_e,fe_e,dfidx_e] = elementType.elementData();

                % Store in the global vectors
                K_ij(counterK+1:counterK+nGleij)     = K_e(:);
                C_ij(counterK+1:counterK+nGleij)     = C_e(:);
                dfidx_ij(counterK+1:counterK+nGleij) = dfidx_e(:);
                fi_i(counterF+1:counterF+nGlei)      = fi_e(:);
                fe_i(counterF+1:counterF+nGlei)      = fe_e(:);
            
                % Update auxiliar variables
                counterK = counterK + nGleij;
                counterF = counterF + nGlei;
                
            end

            % Assemble the matrices and vector
            K      = sparse(iDof,jDof,K_ij);
            C      = sparse(iDof,jDof,C_ij);
            dfidx  = sparse(iDof,jDof,dfidx_ij);
            Fi     = sparse(eDof,ones(this.nDofElemTot,1),fi_i);
            Fe     = sparse(eDof,ones(this.nDofElemTot,1),fe_i);

            % Add contribution of the nodal forces to the external force
            % vector
            Fe = this.addNodalLoad(Fe);
            
        end

        %------------------------------------------------------------------
        % Global system matrices
        function [A,b] = getLinearSystem(this,U,UOld,nonlinearScheme,dt)   

            % Indices for the assembling the matrices and vector
            iDof = zeros(this.sqrNDofElemTot,1);
            jDof = zeros(this.sqrNDofElemTot,1);
            eDof = zeros(this.nDofElemTot,1);

            % Initialize the components of the matrices and vector
            A_ij  = zeros(this.sqrNDofElemTot,1);
            b_i  = zeros(this.nDofElemTot,1);

            % Initialize auxiliar variables
            counterK = 0;
            counterF = 0;

            % Compute and assemble element data
            for el = 1:this.nelem
                % Get the element object
                elementType = this.element(el).type;

                % Update the element displacement vector of each element (TODO: move it to update state variables)
                elementType.DTime = dt;
                elementType.ue    = U(elementType.gle);
                elementType.ueOld = UOld(elementType.gle);

                % Get the vector of the element dof  
                gle_i  = elementType.gle;
                gle_j  = gle_i;
                nGlei  = length(gle_i);
                nGlej  = length(gle_j);
                nGleij = nGlei*nGlej;

                % Get the indices for assemblage
                iDofEl = repmat(gle_i',1,nGlej);
                jDofEl = repmat(gle_j,nGlei,1);
                iDof(counterK+1:counterK+nGleij) = iDofEl(:);
                jDof(counterK+1:counterK+nGleij) = jDofEl(:);
                eDof(counterF+1:counterF+nGlei)  = gle_i';
            
                % Get local matrices and vector
                [A_e,b_e] = elementType.elementLinearSystem(nonlinearScheme);

                % Store in the global vectors
                A_ij(counterK+1:counterK+nGleij) = A_e(:);
                b_i(counterF+1:counterF+nGlei)   = b_e(:);
            
                % Update auxiliar variables
                counterK = counterK + nGleij;
                counterF = counterF + nGlei;
                
            end

            % Assemble the matrices and vector
            A = sparse(iDof,jDof,A_ij);
            b = sparse(eDof,ones(this.nDofElemTot,1),b_i);

            % Add contribution of the nodal forces
            Fe = sparse(this.ndof,1);
            Fe = this.addNodalLoad(Fe);
            b = nonlinearScheme.addNodalForces(b,Fe);

            % Check matrix
            if any(isnan(nonzeros(A)))
                error('Linear system matrix has NaN values');
            end

        end

        %------------------------------------------------------------------
        % Appply the Dirichlet boundary conditions
        function [Aff,bf] = applyDirichletBC(this, A, b, X, nlscheme)
            Aff = A(this.doffree,this.doffree);
            bf  = nlscheme.applyBCtoRHS(A, b, X, this.doffree,this.doffixed);
        end

        %------------------------------------------------------------------
        % Update the state variables from all integration points
        function updateStateVar(this)
            for el = 1:this.nelem
                this.element(el).type.updateStateVar();
            end
        end

        %------------------------------------------------------------------
        % Reorder the nodes to improve computational efficiency
        function resequenceNodes(this)
            % Get auxiliar variables
            nNode   = size(this.NODE,1);
            nElem   = size(this.ELEM,1);
            nNdElem = cellfun(@length,this.ELEM);
            % Size of the connectivity matrix
            nn = sum(nNdElem.^2);
            % Get connectivity matrix
            i=zeros(nn,1); j=zeros(nn,1); s=zeros(nn,1); index=0;
            for el = 1:nElem
              eNode=this.ELEM{el};
              ElemSet=index+1:index+nNdElem(el)^2;
              i(ElemSet) = kron(eNode,ones(nNdElem(el),1))';
              j(ElemSet) = kron(eNode,ones(1,nNdElem(el)))';
              s(ElemSet) = 1;
              index = index + nNdElem(el)^2;
            end
            K = sparse(i,j,s,nNode, nNode);
            % Apply a Symmetric reverse Cuthill-McKee permutation
            p = symrcm(K);
            cNode(p(1:nNode))=1:nNode;
            % Rebuild the nodes and elements matrices
            this.rebuildConnectivity(cNode);
        end

        %------------------------------------------------------------------
        % Updates the connectivity of nodes and elements
        function rebuildConnectivity(this,cNode)
            ELEM_Old = this.ELEM;
            [~,ix,jx] = unique(cNode);
            this.NODE = this.NODE(ix,:);
            for el=1:size(this.ELEM,1)
                for i = 1:length(this.ELEM{el})
                    this.ELEM{el}(i) = jx(ELEM_Old{el}(i));
                end
            end
        end

        %------------------------------------------------------------------
        % Adds pre-existing discontinuities to the model
        function addPreExistingDiscontinuities(this,dSet,additionalData)
            disp("*** Creating the discontinuity elements...");
            if nargin > 2
                this.addDiscontinuityData(additionalData);
            end
            % Check if the mesh is already set
            if (isempty(this.NODE) || isempty(this.ELEM))
                disp('Warning: empty mesh.');
                disp('Warning: the discontinuity set cannot be added.');
                return
            end
            % Create the discontinuity elements
            for i = 1:length(dSet)
                dSet(i).intersectMesh(this) ;
            end
            this.discontinuitySet = dSet;
            this.useEnrichedFormulation(true);
            this.initializeDiscontinuitySegments();
        end

        % -----------------------------------------------------------------
        % Adds some additional discontinuity data
        % To be implemented in the physics whenever required.
        function addDiscontinuityData(~,~)
        end

        %------------------------------------------------------------------
        % Initialize the discontinuity segments
        function initializeDiscontinuitySegments(this)
            if this.useNodalEnrDofs, this.condenseEnrDofs = false; end
            nDiscontinuities = this.getNumberOfDiscontinuities();
            for i = 1:nDiscontinuities
                % Initialize vector with the nodal enrichment dofs
                if this.useNodalEnrDofs
                    nNodesDiscontinuity = size(this.discontinuitySet(i).Xlin,1);
                    nEnrDofs = nNodesDiscontinuity*this.nNodalEnrDofs;
                    nodalEnrDofs = this.ndof+1:(this.ndof+nEnrDofs);
                    this.ndof = this.ndof + nEnrDofs;
                    this.dofenr = [this.dofenr; nodalEnrDofs'];
                    nodalEnrDofs = reshape(nodalEnrDofs, this.nNodalEnrDofs, nNodesDiscontinuity).';
                end
                % Initialize common properties and dofs
                nDiscontinuitySeg = this.discontinuitySet(i).getNumberOfDiscontinuitySegments();
                for j = 1:nDiscontinuitySeg
                    this.discontinuitySet(i).segment(j).t = this.t;
                    if (this.condenseEnrDofs == false) && (this.useNodalEnrDofs == false)
                        this.discontinuitySet(i).segment(j).initializeDofs(this.ndof);
                        this.ndof = this.ndof + this.discontinuitySet(i).segment(j).ndof;
                        this.dofenr = [this.dofenr; this.discontinuitySet(i).segment(j).dof'];
                    end
                    if (this.condenseEnrDofs == false) && (this.useNodalEnrDofs == true)
                        dEnrDofs = [nodalEnrDofs(j,:), nodalEnrDofs(j+1,:)];
                        this.discontinuitySet(i).segment(j).setDofs(dEnrDofs);
                    end
                end
            end
            this.ndofenr = length(this.dofenr);
        end

        %------------------------------------------------------------------
        % Get the total number of discontinuities
        function n = getNumberOfDiscontinuities(this)
            n = length(this.discontinuitySet);
        end

        %------------------------------------------------------------------
        % Flag to use the enriched formulation
        function useEnrichedFormulation(this,flag)
            this.enriched = flag;
        end

        % -----------------------------------------------------------------
        % Print the nodal displacements
        function printResults(this)
            fprintf('\n******** NODAL RESULTS ********\n');
            this.printResultsHeader();
            for i = 1:this.nnodes
                fprintf("  %4d: \t",i);
                for j = 1:this.ndof_nd
                    fprintf("  %+8.4e ",this.U(this.ID(i,j)))
                end
                fprintf("\n");
            end
        end

        %------------------------------------------------------------------
        % Evaluate a field at in a point inside an element
        % Assumes that the point is in the element domain
        function fieldValue = evaluateField(this, field, el, X)
            fieldValue = [];
            if strcmp(field,'Model')
                fieldValue = this.matID(el);
            elseif strcmp(field,'Pressure')
                fieldValue = this.element(el).type.pressureField(X);
            elseif strcmp(field,'Ux')
                u = this.element(el).type.displacementField(X);
                fieldValue = u(1);
            elseif strcmp(field,'Uy')
                u = this.element(el).type.displacementField(X);
                fieldValue = u(2);
            elseif strcmp(field,'E1')
                s = this.element(el).type.strainField(X);
                sp = this.element(el).type.principalStrain(s);
                fieldValue = sp(1);
            elseif strcmp(field,'PEMAG')
                fieldValue = this.element(el).type.plasticstrainMagnitude(X);
            elseif strcmp(field,'Sx')
                s = this.element(el).type.stressField(X);
                fieldValue = s(1);
            elseif strcmp(field,'Sy')
                s = this.element(el).type.stressField(X);
                fieldValue = s(2);
            elseif strcmp(field,'Sxy')
                s = this.element(el).type.stressField(X);
                fieldValue = s(4);
            elseif strcmp(field,'S1')
                s = this.element(el).type.stressField(X);
                sp = this.element(el).type.principalStress(s);
                fieldValue = sp(1);
            elseif strcmp(field,'S2')
                s = this.element(el).type.stressField(X);
                sp = this.element(el).type.principalStress(s);
                fieldValue = sp(2);
            elseif strcmp(field,'Sr')
                s = this.element(el).type.stressField(X);
                sp = this.element(el).type.stressCylindrical(s,X);
                fieldValue = sp(1);
            elseif strcmp(field,'LiquidPressure')
                fieldValue = this.element(el).type.pressureField(X);
            elseif strcmp(field,'CapillaryPressure')
                fieldValue = this.element(el).type.capillaryPressureField(X);
            elseif strcmp(field,'GasPressure')
                fieldValue = this.element(el).type.gasPressureField(X);
            elseif strcmp(field,'LiquidSaturation')
                fieldValue = this.element(el).type.liquidSaturationField(X);
            elseif strcmp(field,'GasSaturation')
                fieldValue = this.element(el).type.gasSaturationField(X);
            end
        end

        %------------------------------------------------------------------
        % Update the result nodes data of each element
        function updateResultVertexData(this,field)
            for el = 1:this.nelem
                this.element(el).type.ue = this.U(this.element(el).type.gle); 
                vertexData = zeros(length(this.element(el).type.result.faces),1);
                for i = 1:length(this.element(el).type.result.faces)
                    X = this.element(el).type.result.vertices(i,:);
                    vertexData(i) = this.evaluateField(field, el, X);
                end
                this.element(el).type.result.setVertexData(vertexData);
            end
        end

        % -----------------------------------------------------------------
        % Plot given field over the mesh
        function plotField(this,field,range,ax)
            if nargin < 3, range = []; end
            if nargin < 4 || isempty(ax)
                figure; 
                ax = gca;
            else
                axes(ax);
                cla(ax);
            end

            this.updateResultVertexData(field)
            FEMPlot(this).plotMesh(ax);
            if strcmp(field,"Model") == false
                if isempty(range)
                    colorbar(ax);
                else
                    clim(ax, range);
                    c = colorbar(ax);
                    c.Limits = range;
                end
            end

        end

        % -----------------------------------------------------------------
        % Plot given field along a given segment
        function plotFieldAlongSegment(this,field, Xi, Xf, npts, axisPlot, ax)
            if nargin < 7 || isempty(ax)
                figure;         % Cria nova figura
                ax = gca;       % Usa o eixo atual
            else
                axes(ax);       % Define o eixo alvo
                cla(ax);        % Limpa o conteúdo
            end
            if nargin < 5
                npts     = 100;
                axisPlot = 'x';
            end
            if nargin < 6
                axisPlot = 'x';
            end
            FEMPlot(this).plotFieldAlongSegment(field, Xi, Xf, npts, axisPlot, ax);
        end

        % -----------------------------------------------------------------
        % Plot field along a given discontinuity
        % Input:
        %   - id: id of the discontinuity 
        function plotFieldAlongDiscontinuiy(this, field, id, axisPlot, ax)
            if nargin < 5 || isempty(ax)
                figure;         
                ax = gca;     
            else
                axes(ax);       
                cla(ax);        
            end
            if nargin < 4
                axisPlot = 'x';
            end

            nDiscontinuitySeg = this.discontinuitySet(id).getNumberOfDiscontinuitySegments();

            % Considering that are two int points per discontinuity segment
            s = zeros(2*nDiscontinuitySeg, 1);
            f = zeros(2*nDiscontinuitySeg, 1);

            % Initial point
            Xi = this.discontinuitySet(id).X(1,:);

            % Fill vectors
            for j = 1:nDiscontinuitySeg
                dof_j = this.discontinuitySet(id).segment(j).dof;
                cElemID = this.discontinuitySet(id).elemID(j);
                cElem = this.element(cElemID).type;
                [Xj, fj] = this.discontinuitySet(id).segment(j).getField(field,this.U(dof_j), cElem);
                DX = Xj - Xi;
                sj = sqrt(DX(:,1).^2 + DX(:,2).^2);
                s(2*j-1:2*j,1) = sj; %Xj(:,2);
                f(2*j-1:2*j,1) = fj;
            end
            
            % Plot
            FEMPlot(this).plotFieldAlongDiscontinuity(s, f, field, axisPlot, ax);
        end

    end
end