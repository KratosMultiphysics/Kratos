%% Shape_LST Class
% This class inherits from the base class 'Shape' to implement the behavior of a linear strain triangular (LST) isoparametric element.
% 
%% Authors
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
%
%% Class Definition
classdef Shape_LST < Shape
    %% Constructor method
    methods
        function this = Shape_LST()
            this = this@Shape('LST');
        end
    end

    %% Public methods
    methods
        %------------------------------------------------------------------
        % Evaluate the shape function at a given point X of a linear 
        % triangular isoparametric element.
         function N = shapeFnc(~,Xn)
            % Natural coordinates of the given point
            r = Xn(1); s = Xn(2);

            % Shape functions
            N1 = 1 - 3*r - 3*s + 4*r*s + 2*r*r + 2*s*s;
            N2 = 2*r*r - r;
            N3 = 2*s*s - s;
            N4 = 4*r - 4*r*r - 4*r*s;
            N5 = 4*r*s;
            N6 = 4*s - 4*r*s - 4*s*s;
            N  = [ N1  N2  N3  N4  N5  N6 ];
         end

         %------------------------------------------------------------------
         % Get the shape function matrix
         function N = shapeFncMtrx(this,Xn)
             % Vector with the shape functions
             N = this.shapeFnc(Xn);
         end

         %------------------------------------------------------------------
         % Get the shape function matrix
         function Nu = NuMtrx(~,N)
             % Vector with the shape functions
             Nu = [N(1)  0.0   N(2)  0.0   N(3)  0.0   N(4)  0.0   N(5)  0.0   N(6)  0.0;
                   0.0   N(1)  0.0   N(2)  0.0   N(3)  0.0   N(4)  0.0   N(5)  0.0   N(6)];
         end

         %------------------------------------------------------------------
         % Compute the derivatives of the shape functions wrt to the
         % natural coordinate s
         function dNdxn = shapeFncDrv(~,Xn)
            % Natural coordinates of the given point
            r = Xn(1); s = Xn(2);

            % Derivatives of the shape functions
            dN1_dr = -3.0 + 4.0*s + 4.0*r;
            dN1_ds = -3.0 + 4.0*r + 4.0*s;
            dN2_dr = 4.0*r - 1.0;
            dN2_ds = 0.0;
            dN3_dr = 0.0;
            dN3_ds = 4.0*s - 1.0;
            dN4_dr = 4.0 - 8.0*r - 4.0*s;
            dN4_ds = -4.0*r;
            dN5_dr = 4.0*s;
            dN5_ds = 4.0*r;
            dN6_dr = -4.0*s;
            dN6_ds = 4.0 - 4.0*r - 8.0*s;

            dNdxn   = [ dN1_dr  dN2_dr  dN3_dr  dN4_dr  dN5_dr  dN6_dr;
                        dN1_ds  dN2_ds  dN3_ds  dN4_ds  dN5_ds  dN6_ds];
         end

         %------------------------------------------------------------------
         % Compute the jacobian matrix
         function J = JacobianMtrx(this,X,Xn)
            % Compute the shape function derivatives wrt to the natural
            % coordinate system
            dNdxn = this.shapeFncDrv(Xn);
            
            % Jacobian matrix
            J = dNdxn*X;
         end

         %------------------------------------------------------------------
         % Compute the determinant of the jacobian
         function detJ = detJacobian(this,X,Xn)
            % Jacobian matrix
            J = this.JacobianMtrx(X,Xn);
            detJ = det(J);
         end

         %------------------------------------------------------------------
         % Compute the derivatives of the shape functions matrix
         function [dNdx,detJ] = dNdxMatrix(this,X,Xn)
            % Jacobian matrix
            J = this.JacobianMtrx(X,Xn);

            % Determinant of the Jacobian matrix
            detJ = det(J);

            % Compute the derivatives of the shape functions wrt to the
            % natural coordinate system
            dNdxn = this.shapeFncDrv(Xn);

            % Compute the derivatives of the shape functions wrt to the
            % global cartesian coordinate system
            dNdx = J\dNdxn;
         end

         %------------------------------------------------------------------
         % Get the integration points:
         % Output:
         %      X: Coordinates of the integration points in the natural
         %         system
         %      w: Weight associated to each integration point
         %      nIntPoints: Total number of integration points
         function [X,w,nIntPoints] = getIntegrationPoints(~,order,~)

             if order == 2
                 X          = [1/6,  1/6;
                               2/3,  1/6;
                               1/6,  2/3]';  
                 w          = [1/6,1/6,1/6];
                 nIntPoints = 3;
             elseif order == 1
                 X          = [1/3, 1/3]';
                 w          = 0.5;
                 nIntPoints = 1;
             end

         end

         %------------------------------------------------------------------
         % Transform a point from the natural coordinate system to the
         % global cartesian coordinate system
         % Input:
         %   NODE : matrix with the x and y coordinates of the nodes of an
         %           CST element
         %   Xn   : vector with the xi and eta coordinates of a point in 
         %          the natural coordinate system
         %
         % Output:
         %   X : vector with the x and y coordinates of a point in the 
         %       global coordinate system
         function X = coordNaturalToCartesian(this,NODE,Xn)

            % Extract the nodal coordinates
            x = NODE(:,1);
            y = NODE(:,2);
            
            % Vector with the shape functions
            Nv = this.shapeFnc(Xn);
            
            % Initialize output
            X = [0.0, 0.0];
            
            % Interpolation the position
            X(1) = Nv(1)*x(1) +  Nv(2)*x(2) +  Nv(3)*x(3) + Nv(4)*x(4) +  Nv(5)*x(5) +  Nv(6)*x(6);
            X(2) = Nv(1)*y(1) +  Nv(2)*y(2) +  Nv(3)*y(3) + Nv(4)*y(4) +  Nv(5)*y(5) +  Nv(6)*y(6);

         end

         % ----------------------------------------------------------------
         % Transform a point from the natural coordinate system to the
         % global cartesian coordinate system
         % Reference: Section 7.3 from the book Concepts and Applications
         % of Finite Element Analysis from Cook et al. (2001)
         function Xn = coordCartesianToNatural(this,NODE,X)
            % Compute the area of the triangle
            A = this.areaTriangle(NODE(1,:), NODE(2,:), NODE(3,:));

            % Compute the areas of the sub-triangles
            A1 = this.areaTriangle(NODE(2,:), NODE(3,:), X);
            A2 = this.areaTriangle(NODE(3,:), NODE(1,:), X);

            % Compute the area coordinates
            xi1 = A1/A;
            xi2 = A2/A;
            xi3 = 1.0 - xi1 - xi2;

            Xn = [xi2, xi3];
         end

        %------------------------------------------------------------------
        % Compute the area of a triangle given the three point in
        % counterclock-wise order
        function A = areaTriangle(~, P1, P2, P3)
            P1P2 = P2 - P1;
            P1P3 = P3 - P1;
            A = (P1P2(1)*P1P3(2) - P1P2(2)*P1P3(1))*0.5;
        end

        %------------------------------------------------------------------
        % Size of the Gram matrix
        % The stress field in a LST element is linear
        function n = dimPolynomialStressInterp(~)
            n = 3;
        end

        %------------------------------------------------------------------
        % Integrand to compute the Gram Matrix
        % The definition of this matrix is associated to the order of the
        % stress field inside the element domain. For a linear
        % triangular isoparametric element, the stress field will be
        % assumed to have a constant
        function dH = integrandGramMtrx(this, node, Xn)
            % Compute the relative coordinate wrt to the centroid of the
            % element
            X    = this.coordNaturalToCartesian(node,Xn);
            Xc   = this.centroid(node);
            Xrel = X - Xc;

            % Gram matrix
            dH = [  1.0          Xrel(1)           Xrel(2);
                   Xrel(1)    Xrel(1)*Xrel(1)   Xrel(2)*Xrel(1);
                   Xrel(2)    Xrel(1)*Xrel(2)   Xrel(2)*Xrel(2)];  
        end

        %------------------------------------------------------------------
        % Size of the stress interpolation vector
        function n = getSizeStressIntVct(~)
            n = 3;
        end

        %------------------------------------------------------------------
        % Size of the stress interpolation vector
        function X = centroid(this,NODE)
            X = coordNaturalToCartesian(this,NODE,[1/3, 1/3]);
        end

        %------------------------------------------------------------------
        % Polynomial stress interpolation
        function p = polynomialStress(~,X)
            p = [1.0; X(1); X(2)];
        end

        %------------------------------------------------------------------
        % Integrand to compute the stress interpolation vector
        function dS = integrandStressIntVct(~,s,Xrel,jumpOrder)
            if jumpOrder == 0
                dS = [  1.0;
                      Xrel(1);
                      Xrel(2)];
            elseif jumpOrder == 1
                dS = [  1.0       s;
                      Xrel(1)  s*Xrel(1);
                      Xrel(2)  s*Xrel(2)];
            end
        end

    end

    methods(Static)
        function [X,w,nIntPoints] = getTriangleLinearQuadrature()
            % Coordinates of the integration points
            X          = [1/6,  1/6;
                          2/3,  1/6;
                          1/6,  2/3]';  
             % Weights
             w          = [1/6,1/6,1/6];
             % Total number of integration points
             nIntPoints = 3;
        end
    end
end