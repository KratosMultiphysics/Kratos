%% RegularElement Class
% This is an abstract class that defines a regular finite element in a finite element mesh.
% It provides properties and methods to define the element's geometry, material properties,
% numerical integration, and other characteristics required for finite element analysis.
%
%% Authors
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
% 
%% Class definition
classdef RegularElement < handle    
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        shape           = [];       % Object of the Shape class
        node            = [];       % Nodes of the fem mesh
        connect         = [];       % Nodes connectivity
        t               = 1.0;      % Thickness
        mat             = [];       % Vector with material properties
        intOrder        = 2;        % Order of the numerical integration
        nnd_el          = 4;        % Number of nodes per element
        ndof_nd         = 1;        % Number of dof per node
        gle             = [];       % Vector of the degrees of freedom
        ngle            = 0;        % Total number of dofs
        ue              = [];       % Element's displacement vector
        ueOld           = [];       % Element's old displacement vector
        due             = [];       % Element's increment displacement
        nIntPoints      = 1;        % Number of integration points
        intPoint        = [];       % Vector with integration point objects
        result          = [];       % Result object to plot the results
        gravityOn       = false;    % Flag to check if the gravity is considered
        g               = 9.806;    % Gravity accelaration (m/s2)
        isEnriched      = false;    % Flag to check if the element is enriched
        massLumping     = false;    % Flag to apply a diagonalization of the compressibility matrix
        lumpStrategy    = 1;        % Id of the diagonalization strategy
        isAxisSymmetric = false;    % Flag to axissymetric models
        DTime           = [];       % Time increment     
        subDivInt       = false;    % Flag to apply a sub-division of the element to define the integration points
    end

    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RegularElement(node,elem,t,mat,intOrder,massLumping,lumpStrategy,isAxisSymmetric)
            if (nargin > 0)
                this.node = node;
                this.nnd_el = size(node,1);
                this.connect = elem;
                this.t = t;
                this.mat = mat;
                this.intOrder = intOrder;
                this.massLumping = massLumping;
                this.lumpStrategy = lumpStrategy;
                this.isAxisSymmetric = isAxisSymmetric;
                order = this.sortCounterClockWise(this.node);
                this.result = Result(this.node(order,:),1:length(this.connect),0.0*ones(this.nnd_el,1),'Model');

                if this.nnd_el == 4
                    this.shape = Shape_ISOQ4();
                elseif this.nnd_el == 8
                    this.shape = Shape_ISOQ8();
                elseif this.nnd_el == 3
                    this.shape = Shape_CST();   
                elseif this.nnd_el == 6
                    this.shape = Shape_LST();
                end
            end
        end
    end

    %% Abstract methods
    methods(Abstract)
        %------------------------------------------------------------------
        % Assemble element matrices and vectors.
        % Outputs:
        %    Ke : element "stiffness" matrix
        %    Ce : element "damping" matrix
        %    fe : element "external force" vector
        %    fi : element "internal force" vector
% dfidu : element matrix of derivative of the internal force with respect to displacement
        [Ke,Ce,fi,fe,dfidu] = elementData(this);
    end

    %% Public methods
    methods
        %------------------------------------------------------------------
        % Assemble element matrices and vectors.
        function [Ae,be] = elementLinearSystem(this,nlscheme)
            [Ke,Ce,fi,fe,dfidu] = this.elementData();
            [Ae,be] = nlscheme.assembleLinearSystem(Ce,Ke,fi,fe,dfidu,this.ue,this.ueOld,this.DTime);
        end
        
        % -----------------------------------------------------------------
        % Update state variables.
        function updateStateVar(this)
            for i = 1:this.nIntPoints
                this.intPoint(i).updateStateVar();
                this.intPoint(i).updateStressVct();
                this.intPoint(i).updateStrainVct();
            end
        end

        %------------------------------------------------------------------
        % Reset state variables.
        function resetIntegrationPts(this)
            for i = 1:this.nIntPoints
                this.intPoint(i).reset();
            end
        end

        %------------------------------------------------------------------
        % Compute element characteristic length.
        function lc = characteristicLength(this)
            lc = this.getDomainArea();
            if strcmp(this.shape.type,'CST') || strcmp(this.shape.type,'LST')
                lc = lc * sqrt(2.0);
            end 
        end

        %------------------------------------------------------------------
        % Compute area of element domain.
        function A = getDomainArea(this)
            A = this.calculateArea(this.node);
        end

        %------------------------------------------------------------------
        % Calculate the area of the element from its vertices.
        function A = calculateArea(~,node)
            % Vertices of the coordinates
            vx = node(:,1); 
            vy = node(:,2);

            % Shifted vertices
            vxS = vx([2:end, 1]);
            vyS = vy([2:end, 1]); 

            % Compute polygon area
            temp = vx.*vyS - vy.*vxS;
            A = 0.5*sum(temp);
        end

        %------------------------------------------------------------------
        % Update result's object vertices property.
        function updateResultVertices(this,configuration)
            if strcmp(configuration,'Deformed')
                Nodes = this.getDeformedConfiguration();
                this.result.setVertices(Nodes);
            end  
        end
    end

    %% Static methods
    methods(Static)
        %------------------------------------------------------------------
        % Sort nodes counterclockwise using the centroid defined by the nodes as reference point.
        function order = sortCounterClockWise(NODE)
            % Centroid coordinate
            cx = mean(NODE(:,1));
            cy = mean(NODE(:,2));

            % Compute the angle that the relative vector of the vertices from the centroid has with the horizontal axis
            a = atan2(NODE(:,2)-cy, NODE(:,1)-cx);

            % Sort the angles
            [~,order] = sort(a);
        end
    end
end
