%% RegularElement_M Class
% This class defines a mechanical finite element. It extends the 
% _RegularElement_ class and provides additional functionality for 
% mechanical analysis, including displacement degrees of freedom, stress 
% and strain computation, and integration point initialization. 
%
%% Methods
% * *initializeIntPoints*: Initializes the integration points for the 
%                          element, including their coordinates, weights, 
%                          and mechanical analysis models.
% * *elementData*: Assembles the element stiffness matrix, damping matrix, 
%                  internal force vector, external force vector, and 
%                  derivative of internal force with respect to 
%                  displacement.
% * *addGravityForces*: Adds the contribution of gravity forces to the 
%                       external force vector.
% * *getNodalDisplacement*: Retrieves the nodal displacement values.
% * *getNodalPressure*: Retrieves the nodal liquid pressure values.
% * *displacementField*: Computes the displacement field at a given 
%                        global Cartesian coordinate.
% * *integrationPointInterpolation*: Computes the interpolation matrix 
%                                    for integration points.
% * *stressField*: Evaluates the stress tensor at a given point by 
%                  extrapolating results from integration points.
% * *strainField*: Evaluates the strain tensor at a given point by 
%                  extrapolating results from integration points.
% * *plasticstrainMagnitude*: Computes the magnitude of the plastic 
%                             strain tensor at a given point.
% * *stressCylindrical*: Transforms the stress tensor to cylindrical 
%                        coordinates.
% * *principalStress*: Computes the principal stresses from the 
%                      stress tensor.
% * *principalStrain*: Computes the principal strains from the 
%                      strain tensor.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
% 
%% Class definition
classdef RegularElement_M < RegularElement    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        glu        = [];            % Displacement dofs
        nglu       = 0;             % Number of regular u-dof
        anm        = 'PlaneStrain'; % Analysis model
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RegularElement_M(node, elem, t, ...
                mat, intOrder, glu, massLumping, lumpStrategy, ...
                isAxisSymmetric,isPlaneStress)
            this = this@RegularElement(node, elem, t, ...
                mat, intOrder, massLumping, lumpStrategy, ...
                isAxisSymmetric);
            this.glu      = glu;
            this.gle      = glu;
            this.nglu     = length(this.glu);
            this.ngle     = length(this.gle);
            if isPlaneStress
                this.anm = 'PlaneStress';
            end
            if isAxisSymmetric
                this.anm = 'AxisSymmetrical';
            end
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Initialize the elements integration points
        function initializeIntPoints(this)

            % Get integration points coordinates and weights
            [X,w,this.nIntPoints] = this.shape.getIntegrationPoints(this.intOrder, this);

            % Get characteristic length
            lc = this.characteristicLength();

            % Initialize the integration points objects
            intPts(this.nIntPoints,1) = IntPoint();
            for i = 1:this.nIntPoints
                constModel = Material_M(this.mat,lc);
                intPts(i) = IntPoint(X(:,i),w(i), constModel);
                intPts(i).initializeMechanicalAnalysisModel(this.anm);
            end
            this.intPoint = intPts;

        end

        %------------------------------------------------------------------
        % This function assembles the element matrices and vectors 
        %
        % Output:
        %    Ke : element "stiffness" matrix
        %    Ce : element "damping" matrix
        %    fe : element "external force" vector
        %    fi : element "internal force" vector
        % dfidu : element matrix of derivative of the internal force with 
        %         respect to displacement
        %
        function [Ke, Ce, fi, fe, dfidu] = elementData(this)

            % Initialize the matrices
            Ke    = zeros(this.nglu, this.nglu);
            Ce    = zeros(this.nglu, this.nglu);
            dfidu = zeros(this.nglu, this.nglu);

            % Initialize external force vector
            fe = zeros(this.nglu, 1);

            % Initialize the internal force vector
            fi = zeros(this.nglu, 1);
            
            % Vector of the nodal dofs
            u  = this.getNodalDisplacement();

            % Numerical integration of the sub-matrices
            for i = 1:this.nIntPoints
                
                % Shape function vector
                N = this.shape.shapeFncMtrx(this.intPoint(i).X);
               
                % Compute the B matrix at the int. point and the detJ
                [dNdx, detJ] = this.shape.dNdxMatrix(this.node,this.intPoint(i).X);

                % Assemble the B-matrix for the mechanical part
                Bu = this.BMatrix(dNdx,N);

                % Compute the strain vector
                this.intPoint(i).strain = Bu * u;

                % Compute the stress vector and the constitutive matrix
                [stress,Duu] = this.intPoint(i).mechanicalLaw();
        
                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;
                if this.isAxisSymmetric
                    c = c * this.shape.axisSymmetricFactor(N,this.node);
                end
                
                % Compute the stiffness sub-matrix
                dfidu = dfidu + Bu' * Duu * Bu * c;

                % Internal force vector
                fi = fi + Bu' * stress * c;
                
                % Compute the gravity forces
                if (this.gravityOn)
                    fe = this.addGravityForces(fe,this.intPoint(i).X,c);
                end
            end
            
        end

        %------------------------------------------------------------------
        % Compute the strain-displacement matrix
        function [B] = BMatrix(this,dNdx,N)
            B = zeros(4,this.nnd_el*2);
            for i = 1:this.nnd_el
                B(1,2*i-1) = dNdx(1,i); 
                B(2,2*i)   = dNdx(2,i);
                B(4,2*i-1) = dNdx(2,i);
                B(4,2*i)   = dNdx(1,i);
                if this.isAxisSymmetric
                    r = N*this.node(:,1);
                    B(3,2*i-1) = N(i)/r;
                end
            end
        end

        %------------------------------------------------------------------
        % Add contribution of the gravity forces to the external force vector
        function fe = addGravityForces(this, fe, Xn, c)

            % Get gravity vector
            grav = this.g * this.mat.porousMedia.b;

            % Shape function matrix
            N  = this.shape.shapeFncMtrx(Xn);
            Nu = this.shape.NuMtrx(N);

            % Get the porous matrix density
            rhos = this.mat.porousMedia.getDensity();

            % Compute the contribution of the gravitational forces
            fe = fe + Nu' * rhos * grav * c;
            
        end

        %------------------------------------------------------------------
        % Compute the forces due to the pore-pressure field
        function fe = porePressureForce(this, pe)

            % Initialize hydro-mechanical coupling matrix. The
            % pore-pressure are discretized at the nodes of the element.
            fe = zeros(this.nglu,1);

            % Identity vector
            m = [1;1;1;0];

            % Numerical integration of fe
            for i = 1:this.nIntPoints

                % Shape function vector
                N = this.shape.shapeFncMtrx(this.intPoint(i).X);

                % Compute the B matrix at the int. point and the detJ
                [dNdx, detJ] = this.shape.dNdxMatrix(this.node,this.intPoint(i).X);

                % Assemble the B-matrix for the mechanical part
                Bu = this.BMatrix(dNdx,N);

                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;
                if this.isAxisSymmetric
                    c = c * this.shape.axisSymmetricFactor(N,this.node);
                end

                % Compute the force due to the pore-pressure
                fe = fe + Bu' * m * N * pe * c;

            end

        end

        %------------------------------------------------------------------
        % Initialize the stresses at the integration point with the given
        % function
        function initialStress(this, stressfnc)
            for i = 1:this.nIntPoints
                % Cartesian coordinates of the integration point
                X = this.shape.coordNaturalToCartesian(this.node,this.intPoint(i).X);
                % Initialize the stresses
                this.intPoint(i).stress    = stressfnc(X);
                this.intPoint(i).stressOld = stressfnc(X);
            end
        end

        %------------------------------------------------------------------
        % Function to get the nodal values of the displacement
        function u = getNodalDisplacement(this)
            u = this.ue(1:this.nglu);
        end

        %------------------------------------------------------------------
        % Function to get old the nodal values of the displacement
        function uOld = getOldNodalDisplacement(this)
            uOld = this.ueOld(1:this.nglu);
        end

        %------------------------------------------------------------------
        % Function to compute the displacement field in the element.
        function u = displacementField(this,X)
        %
        % Input:
        %   X   : position vector in the global cartesian coordinate system
        %
        % Output:
        %   u   : displacement vector evaluated in "X"
        
            % Natural coordinate system
            Xn = this.shape.coordCartesianToNatural(this.node,X);
            
            % Vector with the shape functions
            Nm = this.shape.shapeFncMtrx(Xn);
            Nu = this.shape.NuMtrx(Nm);

            % Displacement dof vector
            uv  = this.getNodalDisplacement();
            
            % Regular displacement field
            u = Nu*uv;
        
        end

        %------------------------------------------------------------------
        % Function to reset the displacements and strains
        function udofs = resetDisplacements(this)

            udofs = this.glu;

            % Reset the displacements
            this.ue(1:this.nglu) = 0.0;
            this.ueOld(1:this.nglu) = 0.0;

            % Reset the strains
            for i = 1:this.nIntPoints
                this.intPoint(i).strain    = zeros(4,1);
                this.intPoint(i).strainOld = zeros(4,1);
            end
        end

        %------------------------------------------------------------------
        % Function to compute the pressure field inside a given element
        function p = pressureField(this,X,pe)
        %
        % Input:
        %   X   : position vector in the global cartesian coordinate system
        %
        % Output:
        %   p   : pressure evaluated in "X"
        
            % Natural coordinate system
            Xn = this.shape.coordCartesianToNatural(this.node,X);
            
            % Vector with the shape functions
            Nm = this.shape.shapeFncMtrx(Xn);

            % capillary field
            p = Nm*pe;
        
        end

        %------------------------------------------------------------------
        % Compute the integration points interpolation matrix
        function S = integrationPointInterpolation(this)
            Q = ones(3,this.nIntPoints);
            for i = 1:this.nIntPoints
                Q(2,i) = this.intPoint(i).X(1);
                Q(3,i) = this.intPoint(i).X(2);
            end
            P = zeros(3,3);
            for i = 1:this.nIntPoints
                P(1,1) = P(1,1) + 1.0;
                P(1,2) = P(1,2) + this.intPoint(i).X(1);
                P(1,3) = P(1,3) + this.intPoint(i).X(2);
                P(2,1) = P(2,1) + this.intPoint(i).X(1);
                P(2,2) = P(2,2) + this.intPoint(i).X(1) * this.intPoint(i).X(1);
                P(2,3) = P(2,3) + this.intPoint(i).X(1) * this.intPoint(i).X(2); 
                P(3,1) = P(3,1) + this.intPoint(i).X(2);
                P(3,2) = P(3,2) + this.intPoint(i).X(2) * this.intPoint(i).X(1);
                P(3,3) = P(3,3) + this.intPoint(i).X(2) * this.intPoint(i).X(2);
            end
            S = P\Q;
        end

        %------------------------------------------------------------------
        % Evaluate the stress tensor in a given point by extrapolating the
        % results from the integration points
        function stress = stressField(this,X,ue)
        %
        % Input:
        %   X   : position vector in the global cartesian coordinate system
            if nargin > 2, this.ue = ue; end
        
            % Natural coordinate system
            Xn = this.shape.coordCartesianToNatural(this.node,X);

            % Get extrapolation matrix
            S = this.integrationPointInterpolation();

            % Matrix with the stress at the integration points
            % Each column corresponds to a stress component:
            % sxx, syy, szz and tauxy
            stressIP = zeros(this.nIntPoints,4);
            for i = 1:this.nIntPoints
                stressIP(i,1) = this.intPoint(i).stress(1);
                stressIP(i,2) = this.intPoint(i).stress(2);
                stressIP(i,3) = this.intPoint(i).stress(3);
                stressIP(i,4) = this.intPoint(i).stress(4);
            end

            % Coefficients for the polynomial approximation 
            c = S * stressIP;

            % Interpolated stress field at the given node
            stress = c' * [1.0 ; Xn(1); Xn(2)];

        end

        %------------------------------------------------------------------
        % Evaluate the strain tensor in a given point by extrapolating the
        % results from the integration points
        function strain = strainField(this,X,ue)
        %
        % Input:
        %   X   : position vector in the global cartesian coordinate system
            if nargin > 2, this.ue = ue; end
        
            % Natural coordinate system
            Xn = this.shape.coordCartesianToNatural(this.node,X);

            % Get extrapolation matrix
            S = this.integrationPointInterpolation();

            % Matrix with the stress at the integration points
            % Each column corresponds to a stress component:
            % sxx, syy, szz and tauxy
            strainIP = zeros(this.nIntPoints,4);
            for i = 1:this.nIntPoints
                strainIP(i,1) = this.intPoint(i).strain(1);
                strainIP(i,2) = this.intPoint(i).strain(2);
                strainIP(i,3) = this.intPoint(i).strain(3);
                strainIP(i,4) = this.intPoint(i).strain(4);
            end

            % Coefficients for the polynomial approximation 
            c = S * strainIP;

            % Interpolated stress field at the given node
            strain = c' * [1.0 ; Xn(1); Xn(2)];

        end

        %------------------------------------------------------------------
        % Evaluate the stress tensor in a given point by extrapolating the
        % results from the integration points
        function pe = plasticstrainMagnitude(this,X,ue)
        %
        % Input:
        %   X   : position vector in the global cartesian coordinate system
            if nargin > 2, this.ue = ue; end
        
            % Natural coordinate system
            Xn = this.shape.coordCartesianToNatural(this.node,X);

            % Get extrapolation matrix
            S = this.integrationPointInterpolation();

            % Matrix with the stress at the integration points
            % Each column corresponds to a stress component:
            % sxx, syy, szz and tauxy
            pstrainIP = zeros(this.nIntPoints,4);
            for i = 1:this.nIntPoints
                pstrainIP(i,1) = this.intPoint(i).plasticstrain(1);
                pstrainIP(i,2) = this.intPoint(i).plasticstrain(2);
                pstrainIP(i,3) = this.intPoint(i).plasticstrain(3);
                pstrainIP(i,4) = this.intPoint(i).plasticstrain(4);
            end

            % Coefficients for the polynomial approximation 
            c = S * pstrainIP;

            % Interpolated stress field at the given node
            pstrain = c' * [1.0 ; Xn(1); Xn(2)];

            pe = norm(pstrain);

        end

        %------------------------------------------------------------------
        % Transforms the stress tensor to cylindrical coordinates
        function sn = stressCylindrical(~,stress,X)

            % Get the stress tensor components
            sx = stress(1);
            sy = stress(2);
            tauxy = stress(4);

            % Compute the angle theta
            theta = atan2(X(2), X(1)); % Angle in radians
            
            % Transform the stresses
            sr = sx * cos(theta)^2 + sy * sin(theta)^2 + 2 * tauxy * cos(theta) * sin(theta);
            stheta = sx * sin(theta)^2 + sy * cos(theta)^2 - 2 * tauxy * cos(theta) * sin(theta);
            taurtheta = (sy - sx) * cos(theta) * sin(theta) + tauxy * (cos(theta)^2 - sin(theta)^2);

            sn = [sr, stheta, taurtheta];
        end

        %------------------------------------------------------------------
        % Function to compute the principal stresses
        function [s1,s2] = principalStress(~,stress)

            % Get the stress tensor components
            sxx = stress(1);
            syy = stress(2);
            sxy = stress(4);

            % Mohr's circle center
            c = (sxx + syy) / 2.0;

            % Mohr's circle radius
            r = sqrt(((sxx - syy)/2.0)^2 + sxy^2);

            % Principal stresses
            s1 = c + r;
            s2 = c - r;

        end

        %------------------------------------------------------------------
        % Function to compute the principal stresses
        function [e1,e2] = principalStrain(~,strain)

            % Get the stress tensor components
            exx = strain(1);
            eyy = strain(2);
            exy = strain(4) / 2.0;

            % Mohr's circle center
            c = (exx + eyy) / 2.0;

            % Mohr's circle radius
            r = sqrt(((exx - eyy)/2.0)^2 + exy^2);

            % Principal stresses
            e1 = c + r;
            e2 = c - r;

        end
    end
end
