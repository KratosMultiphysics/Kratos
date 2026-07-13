%% RegularElement_HM class
% This class defines a single-phase hydromechanical finite element. It 
% extends the _RegularElement_ class and incorporates additional 
% attributes and methods specific to hydromechanical analysis.
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
%                                  matrix based on the element volume and 
%                                  compressibility coefficients.
% * *addGravityForces*: Adds the contribution of gravity forces to the 
%                       external force vector.
% * *getNodalDisplacement*: Retrieves the nodal displacement values.
% * *getNodalPressure*: Retrieves the nodal liquid pressure values.
% * *pressureField*: Computes the pressure field inside the element at a 
%                    given position in the global Cartesian coordinate 
%                    system.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
% 
%% Class definition
classdef RegularElement_HM < RegularElement_M    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        glp        = [];            % Liquid phase pressure dofs
        nglp       = 0;             % Number of regular p-dof
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RegularElement_HM(node, elem, t, ...
                mat, intOrder, glu, glp, massLumping, lumpStrategy, ...
                isAxisSymmetric,isPlaneStress)
            this = this@RegularElement_M(node, elem, t, ...
                mat, intOrder, glu, massLumping, lumpStrategy, ...
                isAxisSymmetric,isPlaneStress);
            this.glp  = glp;
            this.gle  = [glu, glp];
            this.nglp = length(this.glp);
            this.ngle = length(this.gle);
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
                constModel = Material_HM(this.mat);
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

            % Initialize the sub-matrices
            K   = zeros(this.nglu, this.nglu);
            H   = zeros(this.nglp, this.nglp);
            Q   = zeros(this.nglp, this.nglu);
            S   = zeros(this.nglp, this.nglp);

            % Auxiliar zero-matrices and vectors
            Opu = zeros(this.nglp, this.nglu);
            Ouu = zeros(this.nglu, this.nglu);
            Opp = zeros(this.nglp, this.nglp);

            % Initialize external force vector
            feu = zeros(this.nglu, 1);
            fep = zeros(this.nglp, 1);

            % Initialize the internal force vector
            fiu = zeros(this.nglu, 1);
            fip = zeros(this.nglp, 1);
            
            % Vector of the nodal dofs
            u  = this.getNodalDisplacement();
            pl = this.getNodalPressure();

            % Initialize 2D identity vector
            m = [1.0 ; 1.0 ; 1.0 ; 0.0];

            % Initialize the volume of the element
            vol = 0.0;

            % Numerical integration of the sub-matrices
            for i = 1:this.nIntPoints

                % Shape function matrix
                Np = this.shape.shapeFncMtrx(this.intPoint(i).X);
               
                % Compute the B matrix at the int. point and the detJ
                [Bp, detJ] = this.shape.dNdxMatrix(this.node,this.intPoint(i).X);

                % Assemble the B-matrix for the mechanical part
                Bu = this.BMatrix(Bp);

                % Pressure values at the integration point
                pIP = Np * pl;

                % Compute the strain vector
                this.intPoint(i).strain = Bu * u;

                % Compute the stress vector and the constitutive matrix
                [stress,Duu] = this.intPoint(i).mechanicalLaw();
        
                % Compute the permeability matrix
                kh = this.intPoint(i).constitutiveMdl.permeabilityTensor();

                % Get compressibility coefficient
                comp = this.intPoint(i).constitutiveMdl.compressibilityCoeff();

                % Get Biot's coefficient
                biot = this.intPoint(i).constitutiveMdl.biotCoeff();
        
                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;
                if this.isAxisSymmetric
                    c = c * this.shape.axisSymmetricFactor(Np,this.node);
                end
                
                % Compute the stiffness sub-matrix
                K = K + Bu' * Duu * Bu * c;

                % Internal force vector
                fiu = fiu + Bu' * stress * c;

                % Compute the hydromechanical coupling matrices
                Q = Q + Np' * biot * m' * Bu * c;
        
                % Compute permeability sub-matrices
                H = H + Bp' * kh * Bp * c;

                % Compute compressibility matrices
                if ((this.massLumping) && (this.lumpStrategy == 1))
                    S = S + diag(comp*Np*c);
                elseif (this.massLumping == false)
                    S = S + Np' * comp * Np * c;
                end
                
                % Compute the gravity forces
                if (this.gravityOn)
                    [feu,fep] = this.addGravityForces(feu,fep,Np,Bp,kh,pIP,c);
                end

                % Compute the element volume
                vol = vol + c;
            end

            % Compute the lumped mass matrix
            if ((this.massLumping) && (this.lumpStrategy == 2))
                S = lumpedCompressibilityMatrix(this, vol);
            end

            % Assemble the element permeability
            Ke = [ Ouu , -Q';
                   Opu,  H ];

            dfidu = [ K , Opu';
                     Opu, Opp ];
            
            % Assemble the element compressibility matrix
            Ce = [ Ouu , Opu';
                    Q ,   S ];

            % Assemble element internal force vector
            fi = [fiu; fip];

            % Assemble element external force vector
            fe = [feu; fep];
            
        end

        %------------------------------------------------------------------
        % Compute the lumped mass matrices
        function S = lumpedCompressibilityMatrix(this, vol)

            % Get compressibility coefficients
            comp = this.intPoint(1).constitutiveMdl.compressibilityCoeff();

            % Mass distribution factor
            factor = vol / this.nnd_el;

            % Compressibility matrices
            S = comp * factor * eye(this.nglp,this.nglp);

        end

        %------------------------------------------------------------------
        % Add contribution of the gravity forces to the external force vector
        function [feu,fep] = addGravityForces(this,feu,fep,Np,Bp,kh,pl,c)

            % Get gravity vector
            grav = this.g * this.mat.porousMedia.b;

            % Shape function matrix
            Nu = this.shape.NuMtrx(Np);

            % Get the porous matrix density
            rhos = this.mat.porousMedia.getDensity();

            % Get fluid densities
            rhol = this.mat.liquidFluid.getDensity(pl);

            % Compute the contribution of the gravitational forces
            feu = feu + Nu' * rhos * grav * c;
            fep = fep + Bp' * kh * rhol * grav * c;
            
        end

        %------------------------------------------------------------------
        % Function to get the nodal values of the displacement
        function u = getNodalDisplacement(this)
            u = this.ue(1:this.nglu);
        end

        %------------------------------------------------------------------
        % Function to get the nodal values of the liquid pressure
        function pl = getNodalPressure(this)
            a = this.nglu + 1;
            b = this.nglu + this.nglp;
            pl = this.ue(a:b);
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
