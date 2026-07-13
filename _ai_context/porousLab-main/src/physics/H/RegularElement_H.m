%% RegularElement_H Class
% This class defines a single-phase fluid flow finite element. It extends
% the _RegularElement_ class and includes additional properties and 
% methods specific to handling pore pressure degrees of freedom and 
% related computations. 
%
%% Methods
% * *initializeIntPoints*: Initializes the integration points for the 
%                          element using the shape function and material 
%                          properties.
% * *elementData*: Assembles the element stiffness matrix, damping matrix, 
%                  internal force vector, external force vector, and 
%                  derivative of internal force with respect to 
%                  displacement.
% * *lumpedCompressibilityMatrix*: Computes the lumped compressibility 
%                                  matrices based on the element volume 
%                                  and compressibility coefficients.
% * *addGravityForces*: Adds the contribution of gravity forces to the 
%                       external force vector.
% * *getNodalPressure*: Retrieves the nodal values of the liquid pressure.
%
% * *pressureField*: Computes the pressure field at a given position 
%                    inside the element.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
% 
%% Class definition
classdef RegularElement_H < RegularElement    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        glp        = [];            % Pore pressure dofs
        nglp       = 0;             % Number of regular p-dof
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RegularElement_H(node, elem, t, ...
                mat, intOrder, glp, massLumping, lumpStrategy, ...
                isAxisSymmetric)
            this = this@RegularElement(node, elem, t, ...
                mat, intOrder, massLumping, lumpStrategy, ...
                isAxisSymmetric);
            this.glp      = glp;
            this.gle      = glp;
            this.nglp     = length(this.glp);
            this.ngle     = length(this.gle);
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Initialize the elements integration points
        function initializeIntPoints(this)

            % Get integration points coordinates and weights
            [X,w,this.nIntPoints] = this.shape.getIntegrationPoints(this.intOrder);

            % Initialize the integration points objects
            intPts(this.nIntPoints,1) = IntPoint();
            for i = 1:this.nIntPoints
                constModel = Material_H(this.mat);
                intPts(i) = IntPoint(X(:,i),w(i), constModel);
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

            % Initialize the sub-matrices
            Ke    = zeros(this.nglp, this.nglp);
            Ce    = zeros(this.nglp, this.nglp);
            dfidu = zeros(this.nglp, this.nglp);

            % Initialize external force vector
            fe = zeros(this.nglp, 1);
            fi = zeros(this.nglp, 1);
            
            % Vector of the nodal dofs
            pl = this.getNodalPressure();

            % Initialize the volume of the element
            vol = 0.0;

            % Numerical integration of the sub-matrices
            for i = 1:this.nIntPoints

                % Shape function matrix
                Np = this.shape.shapeFncMtrx(this.intPoint(i).X);
               
                % Compute the B matrix at the int. point and the detJ
                [Bp, detJ] = this.shape.dNdxMatrix(this.node,this.intPoint(i).X);

                % Pressure values at the integration point
                pIP = Np * pl;
        
                % Compute the permeability matrix
                kh = this.intPoint(i).constitutiveMdl.permeabilityTensor();

                % Get compressibility coefficient
                comp = this.intPoint(i).constitutiveMdl.compressibilityCoeff();
        
                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;
                if this.isAxisSymmetric
                    c = c * this.shape.axisSymmetricFactor(Np,this.node);
                end
        
                % Compute permeability sub-matrices
                Ke = Ke + Bp' * kh * Bp * c;

                % Compute compressibility matrices
                if ((this.massLumping) && (this.lumpStrategy == 1))
                    Ce = Ce + diag(comp*Np*c);
                elseif (this.massLumping == false)
                    Ce = Ce + Np' * comp * Np * c;
                end
                
                % Compute the gravity forces
                if (this.gravityOn)
                    fe = this.addGravityForces(fe,Bp,kh,pIP,c);
                end

                % Compute the element volume
                vol = vol + c;
            end

            % Compute the lumped mass matrix
            if ((this.massLumping) && (this.lumpStrategy == 2))
                Ce = lumpedCompressibilityMatrix(this, vol);
            end

        end

        %------------------------------------------------------------------
        % Compute the lumped mass matrices
        function S= lumpedCompressibilityMatrix(this, vol)

            % Get compressibility coefficients
            comp = this.intPoint(1).constitutiveMdl.compressibilityCoeff();

            % Mass distribution factor
            factor = vol / this.nnd_el;

            % Compressibility matrices
            S = comp * factor * eye(this.nglp,this.nglp);

        end

        %------------------------------------------------------------------
        % Add contribution of the gravity forces to the external force vector
        function fe = addGravityForces(this,fe,Bp,kh,pl,c)

            % Get gravity vector
            grav = this.g* this.mat.porousMedia.b;

            % Get fluid densities
            rhol = this.mat.liquidFluid.getDensity(pl);

            % Compute the contribution of the gravitational forces
            fe = fe + Bp' * kh * rhol * grav * c;
            
        end

        %------------------------------------------------------------------
        % Function to get the nodal values of the liquid pressure
        function pl = getNodalPressure(this)
            pl = this.ue(1:this.nglp);
        end

        %------------------------------------------------------------------
        % Function to compute the pressure field inside a given element
        function p = pressureField(this,X,ue)
        %
        % Input:
        %   X   : position vector in the global cartesian coordinate system
        %
        % Output:
        %   p   : pressure evaluated in "X"
            if nargin > 2, this.ue = ue; end
        
            % Natural coordinate system
            Xn = this.shape.coordCartesianToNatural(this.node,X);
            
            % Vector with the shape functions
            Nm = this.shape.shapeFncMtrx(Xn);

            % Get nodal pressures
            pl = this.getNodalPressure();

            % capillary field
            p = Nm*pl;
        
        end
    end
end
