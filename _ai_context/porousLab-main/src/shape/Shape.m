%% Shape Class
% This is an abstract class that defines an element shape.
%
%% Authors
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
%
%% Class Definition
classdef Shape < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        type = [];          
    end

    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Shape(type)
            if nargin > 0
                this.type = type;
            end
        end
    end

    %% Abstract methods
    methods (Abstract)
        %------------------------------------------------------------------
         % Evaluate shape function at a given point.
         N = shapeFnc(this,Xn)

         %------------------------------------------------------------------
         % Get shape function matrix.
         N = shapeFncMtrx(this,Xn)

         %------------------------------------------------------------------
         % Get shape function derivatives matrix.
         dNdxi = shapeFncDrv(this,Xn)

         %------------------------------------------------------------------
         % Compute jacobian matrix.
         J = JacobianMtrx(this,X,Xn)

         %------------------------------------------------------------------
         % Compute the determinant of the jacobian.
         detJ = detJacobian(this,X,Xn)

         %------------------------------------------------------------------
         % Transform point coordinates from natural coordinate system to
         % global cartesian coordinate system.
         X = coordNaturalToCartesian(this,NODE,Xn)

         %------------------------------------------------------------------
         % Transform point coordinates from global cartesian coordinate system to
         % natural coordinate system.
         Xn = coordCartesianToNatural(this,NODE,X)

         %------------------------------------------------------------------
         % Get integration points.
         [X,W,n] = getIntegrationPoints(this,intOrder,elem)
    end

    %% Public methods
    methods
        %------------------------------------------------------------------
        % Compute shape centroid.
        function Xc = computeCentroid(~,NODE)
            x = NODE(1:3,1);
            y = NODE(1:3,2);
            polyin = polyshape({x},{y});
            [xc,yc] = centroid(polyin);
            Xc = [xc yc];
        end

        %------------------------------------------------------------------
        % Compute the axisymmetric factor.
        function af = axisSymmetricFactor(~,N,X)
            r = N*X(:,1);
            af = 2.0 * r * pi;
        end
    end
end
