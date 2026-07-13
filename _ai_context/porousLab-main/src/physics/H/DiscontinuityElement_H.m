%% DiscontinuityElement_H Class
% This class defines a hydraulic discontinuity element for modeling 
% discontinuities in porous media mechanics. It inherits from the 
% _DiscontinuityElement_ base class and provides specific implementations 
% for handling hydraulic discontinuities.
%
%% Methods
% * *initializeIntPoints*: Initializes the integration points for the 
%                          discontinuity element. Retrieves integration 
%                          points' coordinates and weights, and creates 
%                          _IntPoint_ objects with the associated material 
%                          model.
% * *elementData*: Computes the element stiffness matrix, internal force 
%                  vector, and other element data based on the input 
%                  displacement vector. Performs numerical integration 
%                  over the integration points.
% 
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class Definition
classdef DiscontinuityElement_H < DiscontinuityElement    
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        ndof_jump = 1;      % Pressure jump dofs
        ndof_int  = 2;      % Discontinuity internal pressure dofs
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = DiscontinuityElement_H(node, mat)
            this = this@DiscontinuityElement(node, mat)
            this.ndof = this.ndof_jump + this.ndof_int;
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        % Initializes the integration points for the element obtaining the
        % coordinates and weights
        function initializeIntPoints(this)

            % Get integration points coordinates and weights
            [X,w,this.nIntPoints] = this.shape.getIntegrationPoints(1);

            % Initialize the integration points objects
            intPts(this.nIntPoints,1) = IntPoint();
            for i = 1:this.nIntPoints
                constModel = MaterialDiscontinuity_H(this.mat);
                intPts(i) = IntPoint(X(:,i),w(i), constModel);
            end
            this.intPoint = intPts;

        end

        %------------------------------------------------------------------
        % Computes the element stiffness matrix, internal force vector and
        % other optional outputs using numerical integration over the
        % element
        % 
        % Outputs:
        %   Ke    - Element stiffness matrix.
        %   Ce    - Element damping matrix.
        %   fi    - Internal force vector.
        %   fe    - External force vector.
        %   dfidu - Derivative of internal force with respect to 
        %           displacement.
        function [Hddi, Sddi, Lcci, Lcji, Lcdi, Ljci, Ljji, Ljdi, Ldci, Ldji, Lddi, fi, fe, dfidu] = elementData(this, celem, id)
            
            % Declare output matrices that won't be used
            fi = []; fe = []; dfidu = [];

            % Initialize the matrices for the numerical integration
            Hddi = zeros(this.ndof_int,this.ndof_int);
            Sddi = zeros(this.ndof_int,this.ndof_int);
            Lcci = zeros(celem.nglp,celem.nglp);
            Lcji = zeros(celem.nglp,this.ndof_jump);
            Lcdi = zeros(celem.nglp,this.ndof_int);
            Ljji = zeros(this.ndof_jump,this.ndof_jump);
            Ljdi = zeros(this.ndof_jump,this.ndof_int);
            Lddi = zeros(this.ndof_int,this.ndof_int);

            % Initialize output matrices
            for i = 1:this.nIntPoints

                % Get the shape function matrix
                N  = this.shape.shapeFnc(this.intPoint(i).X);

                % Cartesian coordinates of the integration point
                Xcar = this.shape.coordNaturalToCartesian(this.node,this.intPoint(i).X);

                % Natural coordinates of the integration in the continuum
                % element
                Xn = celem.shape.coordCartesianToNatural(celem.node,Xcar);

                % Shape function matrix of the continuum
                Np = celem.shape.shapeFncMtrx(Xn);

                % Enriched shape function matrix
                Nenr = celem.enrichedShapeFncValues(id, Np, Xcar);
                Nb = Nenr(2);
                Nt = Nenr(3);

                % Compute the B matrix at the int. point and the detJ
                [dN, detJ] = this.shape.dNdxMatrix(this.node,this.intPoint(i).X);

                % Compute the permeability matrix
                kh = this.intPoint(i).constitutiveMdl.longitudinalPermeability();

                % Get compressibility coefficient
                comp = this.intPoint(i).constitutiveMdl.compressibility();

                % Get leak-offs
                lt = this.intPoint(i).constitutiveMdl.leakoff;
                lb = this.intPoint(i).constitutiveMdl.leakoff;

                % Numerical integration term. The determinant is ld/2.
                c = detJ * this.intPoint(i).w * this.t;

                % Compute the permeability matrix
                Hddi = Hddi + dN' * kh * dN * c;

                % Compute the compressibility matrix
                Sddi = Sddi + N' * comp * N * c;

                % Compute the fluid-flow coupling matrices
                Lcci = Lcci + (lt + lb) * (Np' * Np) * c;
                Lcji = Lcji + Np' * (lt * Nt + lb * Nb) * c;
                Lcdi = Lcdi + (lt + lb) * Np' * N * c;
                Ljji = Ljji + (lt * Nt * Nt + lb * Nb * Nb) * c;
                Ljdi = Ljdi + (lt * Nt * N + lb * Nb * N) * c;
                Lddi = Lddi + (lt + lb) * (N' * N) * c;
            end
            Ljci = Lcji';
            Ldci = Lcdi';
            Ldji = Ljdi';
        end
    end
end