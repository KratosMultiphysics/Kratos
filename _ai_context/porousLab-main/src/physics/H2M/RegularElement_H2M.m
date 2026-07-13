%% RegularElement_H2M Class
% This class defines a finite element for a two-phase flow formulation 
% using displacements (ux, uy), liquid pressure (Pl) and gas pressure (Pg) 
% as primary variables. It extends the _RegularElement_ class and 
% incorporates additional attributes and methods specific to 
% hydromechanical coupling in porous media.
%
%% Methods
% * *initializeIntPoints*: Initializes the integration points for the 
%                          element using the shape function and material 
%                          properties.
% * *elementData*: Assembles the coupled mechanical, liquid-pressure, and
%                  gas-pressure residual and tangent terms, including
%                  deformation, advective flow, gravity, storage,
%                  saturation, density, and relative permeability
%                  contributions.
% * *getNodalDisplacement*: Retrieves the nodal displacement values.
% * *getNodalLiquidPressure*: Retrieves the nodal liquid pressure values.
% * *getNodalGasPressure*: Retrieves the nodal gas pressure values.
% * *getNodalCapillaryPressure*: Retrieves the nodal capillary pressure 
%                                values.
% * *pressureField*: Computes the pressure field at a given position 
%                    inside the element.
% * *gasPressureField*: Computes the gas pressure field at a given 
%                       position inside the element.
% * *capillaryPressureField*: Computes the capillary pressure field at a 
%                             given position inside the element.
% * *liquidSaturationField*: Computes the liquid saturation field at a 
%                            given position inside the element.
% * *gasSaturationField*: Computes the gas saturation field at a given 
%                         position inside the element.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
% 
%% Class definition
classdef RegularElement_H2M < RegularElement_M    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        glp        = [];            % Liquid phase pressure dofs
        glpg       = [];            % Gas phase pressure dofs
        nglp       = 0;             % Number of regular p-dof
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RegularElement_H2M(node, elem, t, ...
                mat, intOrder, glu, glp, glpg, massLumping, lumpStrategy, ...
                isAxisSymmetric,isPlaneStress)
            this = this@RegularElement_M(node, elem, t, ...
                mat, intOrder, glu, massLumping, lumpStrategy, ...
                isAxisSymmetric,isPlaneStress);
            this.glp      = glp;
            this.glpg     = glpg;
            this.gle      = [glu, glp, glpg];
            if (length(this.glp) ~= length(this.glpg))
                error('Wrong number of pressure dofs');
            end
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
                constModel = Material_H2M(this.mat);
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

            % Get constitutive model
            constModel = this.intPoint(1).constitutiveMdl;

            % Get gravity vector
            grav = this.g * this.mat.porousMedia.b;

            % Get the porous matrix density
            rhos = this.mat.porousMedia.getDensity();

            % Get the fluid viscosity
            mul = this.mat.liquidFluid.mu;
            mug = this.mat.gasFluid.mu;

            % Get porosity
            phi = constModel.porousMedia.phi;

            % Get the fluids bulk modulus
            Klb  = constModel.liquidFluid.getBulkModulus();
            Kgb  = 1e25; %constModel.gasFluid.getBulkModulus();
            
            % Vector of the nodal dofs
            u  = this.getNodalDisplacement();
            pl = this.getNodalLiquidPressure();
            pg = this.getNodalGasPressure();
            pc = pg - pl;

            % Vector with the old nodal dofs
            plOld = this.getOldNodalLiquidPressure(); 
            pgOld = this.getOldNodalGasPressure();
            pcOld = pgOld - plOld;

            % Fill nodal state variables
            Sl      = zeros(this.nnd_el,1);
            SlOld   = zeros(this.nnd_el,1);
            dSldpc  = zeros(this.nnd_el,1);
            rhol    = zeros(this.nnd_el,1);
            rhog    = zeros(this.nnd_el,1);
            rholOld = zeros(this.nnd_el,1);
            rhogOld = zeros(this.nnd_el,1);
            for i = 1:this.nnd_el
                % Liquid saturation degree
                Sl(i)      = constModel.saturationDegree(pc(i));
                SlOld(i)   = constModel.saturationDegree(pcOld(i));
                % Liquid saturation derivative wrt pc
                dSldpc(i)  = constModel.derivativeSaturationDegree(pc(i));
                % Get fluid densities
                rhol(i)    = this.mat.liquidFluid.getDensity(pl(i));
                rholOld(i) = this.mat.liquidFluid.getDensity(plOld(i));
                rhog(i)    = this.mat.gasFluid.getDensity(pg(i));
                rhogOld(i) = this.mat.gasFluid.getDensity(pgOld(i));
            end

            % Auxiliary derivatives
            dpcdpg = 1.0;
            dpcdpl = -1.0;
            
            % Derivative of the mean saturation wrt the nodal saturation
            dSlmSli = 1.0/this.nnd_el;

            % Gas saturation
            Sg = 1.0 - Sl;
            SgOld = 1.0 - SlOld;
            
            % Derivative of the mean liquid saturation wrt the pressure
            dSlmdpl = dSlmSli * dSldpc * dpcdpl;
            dSlmdpg = dSlmSli * dSldpc * dpcdpg;

            % Nodal mass increment
            Dml = Sl .* rhol - SlOld .* rholOld;
            Dmg = Sg .* rhog - SgOld .* rhogOld;

            % Compute the relative permeability
            [klr, kgr] = constModel.relativePermeabilities(mean(Sl));

            % Derivative of the relative permeability wrt the saturation
            [dklrdSlm, dkgrdSlm] = constModel.derivativeRelPerm(mean(Sl));

            % Derivative of the relative permeability wrt the pressure
            dklrdPl = dklrdSlm * dSlmdpl;
            dklrdPg = dklrdSlm * dSlmdpg;
            dkgrdPl = dkgrdSlm * dSlmdpl;
            dkgrdPg = dkgrdSlm * dSlmdpg;

            % Average bulk density
            rhoavg = (1.0 - phi) * rhos + phi * (mean(Sl) * mean(rhol) + mean(Sg) * mean(rhog));

            % Initialize 2D identity vector
            m = [1.0 ; 1.0 ; 1.0 ; 0.0];

            % Initialize the volume of the element
            vol = 0.0;

            % Numerical integration terms
            fiu = zeros(this.nglu,1);
            dfiudu = zeros(this.nglu,this.nglu);
            H = zeros(this.nglp, this.nglp);
            Q = zeros(this.nglu, this.nglp);
            fgravu = zeros(this.nglu, 1);
            fgravp = zeros(this.nglp, 1);
            fv = zeros(this.nglp, 1);
            for i = 1:this.nIntPoints

                % Shape function matrix
                Np = this.shape.shapeFncMtrx(this.intPoint(i).X);

                % Shape function matrix to interpolate the displacements
                Nu = this.shape.NuMtrx(Np);
               
                % Compute the B matrix at the int. point and the detJ
                [Bp, detJ] = this.shape.dNdxMatrix(this.node,this.intPoint(i).X);

                % Assemble the B-matrix for the mechanical part
                Bu = this.BMatrix(Bp);
               
                % Compute the strain vector
                this.intPoint(i).strain = Bu * u;

                % Compute the stress vector and the constitutive matrix
                [stress,Duu] = this.intPoint(i).mechanicalLaw();

                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;
                if this.isAxisSymmetric
                    c = c * this.shape.axisSymmetricFactor(Np,this.node);
                end
                
                % Compute the stiffness sub-matrix
                dfiudu = dfiudu + Bu' * Duu * Bu * c;

                % Internal force vector
                fiu = fiu + Bu' * stress * c;

                % Get porous media and fluid parameters
                K = this.mat.porousMedia.intrinsicPermeabilityMatrix();

                % Compute permeability matrix
                H = H + Bp' * K * Bp * c;
                
                % Gravity force
                fgravu = fgravu + Nu' * rhoavg * grav * c;
                fgravp = fgravp + Bp' * K * grav * c;

                % Compute the hydromechanical coupling matrix
                Q = Q + Bu' * m * Np * c;
                
                % Volumetric strain rate
                devdt = m' * (this.intPoint(i).strain - this.intPoint(i).strainOld) / this.DTime;

                % Volumetric strain forces
                fv = fv + Np' * devdt * c;

                % Compute the element volume
                vol = vol + c;
            end

            % Mechanical internal force: add contribution of the
            % pore-pressures
            fiu = fiu - mean(Sl) * Q * pl - mean(Sg) * Q * pg;

            % Derivatives of the mechanical internal force wrt to the
            % pressures
            dfiudpl = Q * pc * dSlmdpl' - mean(Sl) * Q;
            dfiudpg = Q * pc * dSlmdpg' - mean(Sg) * Q;

            % Advective forces
            fil = (klr / mul) * H * pl;
            fig = (kgr / mug) * H * pg;

            % Derivatives of the advective forces
            Hll = (klr / mul) * H;
            Hgg = (kgr / mug) * H;
            Hll = Hll + (H * pl) * dklrdPl' /mul;
            Hlg = (H * pl) * dklrdPg' /mul;
            Hgl = (H * pg) * dkgrdPl' /mug;
            Hgg = Hgg + (H * pg) * dkgrdPg' /mug;

            % Volumetric strain forces
            fil = fil + mean(Sl) * fv;
            fig = fig + mean(Sg) * fv;

            % Derivatives of the volumetric strain forces
            dfildu = mean(Sl) * Q' / this.DTime;
            dfigdu = mean(Sg) * Q' / this.DTime;
            Hll = Hll + fv * dSlmdpl';
            Hlg = Hlg + fv * dSlmdpg';
            Hgl = Hgl - fv * dSlmdpl';
            Hgg = Hgg - fv * dSlmdpg';

            % Storage terms
            masscoeff = phi * (vol / this.nnd_el) / this.DTime;
            fil = fil + (Dml ./ rhol) * masscoeff;
            fig = fig + (Dmg ./ rhog) * masscoeff;
            Cll = ((SlOld .* rholOld ./ rhol) / Klb - dSldpc) * masscoeff;
            Clg = (dSldpc) * masscoeff;
            Cgg = ((SgOld .* rhogOld ./ rhog) / Kgb - dSldpc) * masscoeff;

            % Derivatives of the residual of the fluid flow equations
            dfildpl = Hll + diag(Cll);
            dfildpg = Hlg + diag(Clg);
            dfigdpl = Hgl + diag(Clg);
            dfigdpg = Hgg + diag(Cgg);

            % Compute the gravity forces
            if (this.gravityOn)
                fiu = fiu - fgravu;
                fil = fil - (klr / mul) * mean(rhol) * fgravp;
                fig = fig - (kgr / mug) * mean(rhog) * fgravp;
                dfildpl = dfildpl - (1.0/mul) * fgravp * (mean(rhol) * dklrdPl + (klr * dSlmSli / Klb) * rhol)';
                dfildpg = dfildpg - (1.0/mul) * fgravp * (mean(rhol) * dklrdPg)';
                dfigdpl = dfigdpl - (1.0/mug) * fgravp * (mean(rhog) * dkgrdPl)';
                dfigdpg = dfigdpg - (1.0/mug) * fgravp * (mean(rhog) * dkgrdPg + (kgr * dSlmSli / Kgb) * rhog)';
            end
            
            % Assemble the element matrices
            dfidu = [ dfiudu, dfiudpl , dfiudpg;
                      dfildu, dfildpl , dfildpg;
                      dfigdu, dfigdpl , dfigdpg ];

            % Assemble element internal force vector
            fi = [fiu; fil; fig];

            % Initialize matrices that are not being used
            fe = zeros(this.ngle, 1);
            Ce = zeros(this.ngle, this.ngle);
            Ke = zeros(this.ngle, this.ngle);
            
        end

        %------------------------------------------------------------------
        % Function to get the nodal values of the liquid pressure
        function pl = getNodalLiquidPressure(this)
            a = this.nglu + 1;
            b = this.nglu + this.nglp;
            pl = this.ue(a:b);
        end

        %------------------------------------------------------------------
        % Function to get the old nodal values of the liquid pressure
        function plOld = getOldNodalLiquidPressure(this)
            a = this.nglu + 1;
            b = this.nglu + this.nglp;
            plOld = this.ueOld(a:b);
        end

        %------------------------------------------------------------------
        % Function to get the nodal values of the gas pressure
        function pg = getNodalGasPressure(this)
            a = this.nglu + this.nglp + 1;
            pg = this.ue(a:end);
        end

        %------------------------------------------------------------------
        % Function to get the old nodal values of the gas pressure
        function pgOld = getOldNodalGasPressure(this)
            a = this.nglu + this.nglp + 1;
            pgOld = this.ueOld(a:end);
        end

        %------------------------------------------------------------------
        % Function to get the nodal values of the capillary pressure
        function pc = getNodalCapillaryPressure(this)
            pl = this.getNodalLiquidPressure();
            pg = this.getNodalGasPressure();
            pc = pg - pl;
        end

        %------------------------------------------------------------------
        % Function to get the old nodal values of the capillary pressure
        function pcOld = getOldNodalCapillaryPressure(this)
            plOld = this.getOldNodalLiquidPressure();
            pgOld = this.getOldNodalGasPressure();
            pcOld = pgOld - plOld;
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
            pl = this.getNodalLiquidPressure();

            % capillary field
            p = Nm*pl;
        
        end

        %------------------------------------------------------------------
        % Function to compute the pressure field inside a given element
        function p = gasPressureField(this,X,ue)
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
            pg = this.getNodalGasPressure();

            % capillary field
            p = Nm*pg;
        
        end

        %------------------------------------------------------------------
        % Function to compute the pressure field inside a given element
        function p = capillaryPressureField(this,X,ue)
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
            pc = this.getNodalCapillaryPressure();

            % capillary field
            p = Nm*pc;
        
        end

        %------------------------------------------------------------------
        % Function to compute the liquid saturation field inside a given element
        function Sl = liquidSaturationField(this,X,ue)

            if nargin > 2, this.ue = ue; end
        
            % Natural coordinate system
            Xn = this.shape.coordCartesianToNatural(this.node,X);
            
            % Vector with the shape functions
            Nm = this.shape.shapeFncMtrx(Xn);

            % Capillary pressure at the given point
            pc = Nm*this.getNodalCapillaryPressure();

            % Compute the liquid saturation degree
            Sl = this.intPoint(1).constitutiveMdl.saturationDegree(pc);
        
        end

        %------------------------------------------------------------------
        % Function to compute the gas saturation field inside a given element
        function Sg = gasSaturationField(this,X,ue)
            
            if nargin > 2, this.ue = ue; end
        
            % Natural coordinate system
            Xn = this.shape.coordCartesianToNatural(this.node,X);
            
            % Vector with the shape functions
            Nm = this.shape.shapeFncMtrx(Xn);

            % Capillary pressure at the given point
            pc = Nm*this.getNodalCapillaryPressure();

            % Compute the liquid saturation degree
            Sl = this.intPoint(1).constitutiveMdl.saturationDegree(pc);

            % Gas saturation degree
            Sg = 1.0 - Sl;
        
        end
    end
end
