%% DiscontinuityElement_M Class
% This class defines a mechanical discontinuity element that extends the
% _DiscontinuityElement_ base class. It includes additional functionality
% for handling stretching and relative rotation modes, as well as methods
% for initializing integration points and computing element data.
%
%% Methods
% * *addStretchingMode*: Enables or disables the stretching mode based on 
%                        the flag.
% * *addRelRotationMode*: Enables or disables the relative rotation mode 
%                         based on the flag.
% * *initializeIntPoints*: Initializes the integration points for the 
%                          discontinuity element. Retrieves integration 
%                          points' coordinates and weights, and creates 
%                          _IntPoint_ objects with the associated material 
%                          model.
% * *elementData*: Computes the element stiffness matrix, internal force 
%                  vector, and other element data based on the input 
%                  displacement vector. Performs numerical integration 
%                  over the integration points.
% * *enrichmentInterpolationMatrix*: Computes the enrichment interpolation 
%                                    matrix for the given integration point
%                                    coordinates, reference point, and 
%                                    tangential vector.
% 
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00: Initial version (January 2024).
%
%% Class Definition
classdef DiscontinuityElement_M < DiscontinuityElement    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        tangentialStretchingMode = false;
        normalStretchingMode     = false;
        relRotationMode          = false;
        DP                       = [];
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = DiscontinuityElement_M(node, mat)
            this = this@DiscontinuityElement(node, mat)
            this.ndof = 2;
            this.nNodalDofs = 4;
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        % Enables the stretching mode. If enables, the number of degrees of
        % freedom increases by 1 or 2
        function addStretchingMode(this,flagTangential, flagNormal)
            this.tangentialStretchingMode = flagTangential;
            this.normalStretchingMode = flagNormal;
            if flagTangential == true
                this.ndof = this.ndof + 1;
            end
            if flagNormal == true
                this.ndof = this.ndof + 1;
            end
        end

        %------------------------------------------------------------------
        % Enables the rotation mode. If enables, the number of degrees of
        % freedom increases by 1
        function addRelRotationMode(this,flag)
            this.relRotationMode = flag;
            if flag == true
                this.ndof = this.ndof + 1;
            end
        end

        %------------------------------------------------------------------
        % Function to reset the displacements and strains
        function resetDisplacements(this)
            % Reset the strains
            for i = 1:this.nIntPoints
                this.intPoint(i).strain    = zeros(2,1);
                this.intPoint(i).strainOld = zeros(2,1);
            end
        end

        %------------------------------------------------------------------
        % Displacement jump order
        function n = displacementJumpOrder(this)
            n = 0;
            if (this.tangentialStretchingMode || this.normalStretchingMode || this.relRotationMode) 
                n = 1;
            end
        end

        %------------------------------------------------------------------
        % Matrix that transforms the nodal dofs into the local dofs
        function MR = getDofTransformationMtrx(this)
            if this.useNodalEnrDofs == false
                MR = eye(this.ndof);
            else
                ld = this.ld();
                M = [ 0.5    ,  0.0    , 0.5    , 0.0;
                      0.0    ,  0.5    , 0.0    , 0.5;
                     -1.0/ld ,  0.0    , 1.0/ld , 0.0;
                      0.0    , -1.0/ld , 0.0    , 1.0/ld];
                mn = this.rotationFromGlobalToLocal();
                R  = blkdiag(mn,mn);
                MR = M * R;
                lin = [1,2];
                if this.tangentialStretchingMode
                    lin = [lin,3];
                end
                if this.relRotationMode
                    lin = [lin,4];
                end
                MR = MR(lin,:);
            end
        end
        
        %------------------------------------------------------------------
        % Initializes the integration points for the element obtaining the
        % coordinates and weights
        function initializeIntPoints(this)

            % Get integration points coordinates and weights
            [X,w,this.nIntPoints] = this.shape.getIntegrationPoints(1);

            % Initialize the integration points objects
            intPts(this.nIntPoints,1) = IntPoint();
            for i = 1:this.nIntPoints
                constModel = MaterialDiscontinuity_M(this.mat);
                intPts(i) = IntPoint(X(:,i),w(i), constModel);
                intPts(i).initializeMechanicalAnalysisModel('Interface');
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
        function [Ke, Ce, fi, fe, dfidu] = elementData(this, ae)
            
            % Declare output matrices that won't be used
            Ce = []; fe = []; dfidu = [];

            % Initialize the matrices for the numerical integration
            if this.useNodalEnrDofs
                ndofs = this.nNodalDofs;
            else
                ndofs = this.ndof;
            end
            Ke = zeros(ndofs,ndofs);
            fi = zeros(ndofs,1);

            % Get the lenght of the discontinuity
            ld = this.ld();

            % Get the discontinuity reference point
            Xr = this.referencePoint();

            % Get the discontinuity tangential vector
            m = this.tangentialVector();

            % Transformation matrix to use the nodal degrees of freedom
            M = this.getDofTransformationMtrx();

            % Initialize output matrices
            for i = 1:this.nIntPoints

                % Cartesian coordinates of the integration point 
                X = this.shape.coordNaturalToCartesian(this.node,this.intPoint(i).X);

                % Get the shape function matrix
                Ndl = this.enrichmentInterpolationMatrix(X,Xr,m);
                if this.useNodalEnrDofs
                    Nd = this.nodalEnrichmentInterpolationMatrix(this.intPoint(i).X);
                else
                    Nd = Ndl;
                end

                % Compute the strain vector
                this.intPoint(i).strain = Ndl * M * ae;

                % Compute the stress vector and the constitutive matrix
                [td,Td] = this.intPoint(i).mechanicalLaw();

                % Numerical integration term. The determinant is ld/2.
                c = 0.5 * ld * this.intPoint(i).w * this.t;

                % Compute the stiffness sub-matrix
                Ke = Ke + (Nd' * Td * Nd) * c;

                % Compute the internal force vector
                fi = fi + (Nd' * td) * c;

            end
        end

        %------------------------------------------------------------------
        % Computes the enrichment interpolation matrix for the given
        % integration point coordinates, reference point and tangential
        % vector
        function Nd = enrichmentInterpolationMatrix(this,X,Xr,m)
            Nd = zeros(2,this.ndof);
            Nd(1,1) = 1.0;
            Nd(2,2) = 1.0;
            if this.ndof > 2
                s = m' * (X' - Xr');
                c = 3;
                if this.tangentialStretchingMode
                    Nd(1,c) = s;
                    c = c + 1;
                end
                if this.relRotationMode
                    Nd(2,c) = s;
                end
            end
        end

        %------------------------------------------------------------------
        % Computes the enrichment interpolation matrix for the given
        % integration point coordinate. This shape function is used
        % whenever nodal dofs are used. The direct use of this shape
        % function avoid undesired coupling terms in the dofs that appear 
        % from the transformation from global to local dofs. 
        function Nd = nodalEnrichmentInterpolationMatrix(this,xi)
            mn = this.rotationFromGlobalToLocal();
            R  = blkdiag(mn,mn);
            Nd = 0.5*[1.0-xi , 0.0   , 1.0+xi, 0.0;
                      0.0    , 1.0-xi, 0.0   , 1.0+xi]*R;
        end

        %------------------------------------------------------------------
        % Projection matrix
        function P = projectionMatrix(this)
            n = this.normalVector();
            P = [n(1) , 0.0;
                 0.0  , n(2);
                 0.0  , 0.0;
                 n(2) , n(1)];
        end

        %------------------------------------------------------------------
        % Integrate the polynomial stress interpolation of the continuum
        % along the discontinuity
        function S = intPolynomialStressIntp(this, celem)

            % Initialize variables
            jumpOrder     = this.displacementJumpOrder();
            dimPolyStress = celem.shape.dimPolynomialStressInterp();
            S = zeros(dimPolyStress, jumpOrder + 1);

            % Get the discontinuity geometric properties
            Xr = this.referencePoint();
            m  = this.tangentialVector();
            ld = this.ld();

            % Numerical integration
            for i = 1:this.nIntPoints

                % Numerical integration term. The determinant is ld/2.
                c = 0.5 * ld * this.intPoint(i).w * this.t;

                % Cartesian coordinates of the integration point 
                X = this.shape.coordNaturalToCartesian(this.node,this.intPoint(i).X);

                % Continuum stress interpolation polynomial
                p = celem.shape.polynomialStress(X);

                if (jumpOrder == 0)
                    S = S + p * c;
                elseif (jumpOrder == 1)
                    % Tangential coordinate
                    s = m' * (X' - Xr');
                    S = S + [p , s*p] * c;
                end
            end
        end

        %------------------------------------------------------------------
        % Get specified field. Fill the coordinate matrix and the field.
        function [X, f] = getField(this,field,~,~)
            if strcmp(field,'Sn')
                [X, f] = this.getCohesiveStresses(2); 
            elseif strcmp(field,'St')
                [X, f] = this.getCohesiveStresses(1);
            elseif strcmp(field,'Dn')
                [X, f] = this.getDisplacementJump(2);
            elseif strcmp(field,'Dt')
                [X, f] = this.getDisplacementJump(1);
            elseif strcmp(field,'Aperture')
                [X, f] = this.getAperture();
            else 
                X = []; f = [];
            end
        end

        %------------------------------------------------------------------
        % Get cohesive stresses component
        function [X, f] = getCohesiveStresses(this,stressId)
            X = zeros(this.nIntPoints,2);
            f = zeros(this.nIntPoints,1);
            for i = 1:this.nIntPoints
                % Cartesian coordinates of the integration point 
                X(i,:) = this.shape.coordNaturalToCartesian(this.node,this.intPoint(i).X);
                % Get cohesive stresses
                f(i) = this.intPoint(i).stress(stressId);
            end
        end

        %------------------------------------------------------------------
        % Get cohesive stresses component
        function [X, f] = getAperture(this)
            X = zeros(this.nIntPoints,2);
            f = zeros(this.nIntPoints,1);
            for i = 1:this.nIntPoints
                % Cartesian coordinates of the integration point 
                X(i,:) = this.shape.coordNaturalToCartesian(this.node,this.intPoint(i).X);
                % Get aperture component
                f(i) = 1.0E3 * this.intPoint(i).statevar(1); % Convert to mm
            end
        end

        %------------------------------------------------------------------
        % Get cohesive stresses component
        function [X, f] = getDisplacementJump(this,strainId)
            X = zeros(this.nIntPoints,2);
            f = zeros(this.nIntPoints,1);
            for i = 1:this.nIntPoints
                % Cartesian coordinates of the integration point 
                X(i,:) = this.shape.coordNaturalToCartesian(this.node,this.intPoint(i).X);
                % Get strain component
                f(i) = this.intPoint(i).strain(strainId);
            end
        end

        %------------------------------------------------------------------
        function fe = porePressureForce(this, pd)

            % Initialize the matrices for the numerical integration
            fe = zeros(this.ndof,1);

            % Get the lenght of the discontinuity
            ld = this.ld();

            % Get the discontinuity geometry
            Xr = this.referencePoint();
            m  = this.tangentialVector();
            n  = [0;1];

            % Initialize output matrices
            for i = 1:this.nIntPoints

                % Shape function vector
                N = this.shape.shapeFnc(this.intPoint(i).X);

                % Cartesian coordinates of the integration point 
                X = this.shape.coordNaturalToCartesian(this.node,this.intPoint(i).X);

                % Get the shape function matrix
                Nd = this.enrichmentInterpolationMatrix(X,Xr,m);

                % Numerical integration term. The determinant is ld/2.
                c = 0.5 * ld * this.intPoint(i).w * this.t;

                % Pore-pressure at the integration point
                pdi = N * pd;

                % Compute the internal force vector
                fe = fe + Nd' * n * pdi * c;

            end
        end

    end
end