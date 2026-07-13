%% RegularElement_H2 Class
% This class defines a finite element for a two-phase flow formulation 
% using the liquid pressure (Pl) and the gas pressure (Pg) as primary 
% variables. It extends the _RegularElement_ class and provides methods 
% for initializing integration points, assembling element matrices and 
% vectors, and computing various fields such as pressure and saturation.
%
%% Methods
% * *initializeIntPoints*: Initializes the integration points for the 
%                          element using the shape function and material 
%                          properties.
% * *elementData*: Assembles the element residual and tangent terms for
%                  liquid pressure and gas pressure, including advective,
%                  gravity, storage, saturation, density, and relative
%                  permeability contributions.
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
classdef RegularElement_H2 < RegularElement    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        glp        = [];            
        glpg       = [];            % Vector of the regular degrees of freedom
        nglp       = 0;             % Number of regular p-dof
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RegularElement_H2(node, elem, t, ...
                mat, intOrder, glp, glpg, massLumping, lumpStrategy, ...
                isAxisSymmetric)
            this = this@RegularElement(node, elem, t, ...
                mat, intOrder, massLumping, lumpStrategy, ...
                isAxisSymmetric);
            this.glp      = glp;
            this.glpg     = glpg;
            this.gle      = [glp , glpg];
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
                constModel = Material_H2(this.mat);
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

            % Get constitutive model
            constModel = this.intPoint(1).constitutiveMdl;

            % Get gravity vector
            grav = this.g * this.mat.porousMedia.b;

            % Get the fluid viscosity
            mul = this.mat.liquidFluid.mu;
            mug = this.mat.gasFluid.mu;

            % Get porosity
            phi = constModel.porousMedia.phi;

            % Get the fluids bulk modulus
            Klb  = constModel.liquidFluid.K;
            Kgb  = constModel.gasFluid.K;
            
            % Vector of the nodal pore-pressure dofs
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
            
            % Derivative of the mean saturation wrt the nodal saturation
            dSlmSli = 1.0/this.nnd_el;

            % Gas saturation
            Sg = 1.0 - Sl;
            SgOld = 1.0 - SlOld;

            % Nodal mass increment
            Dml = Sl .* rhol - SlOld .* rholOld;
            Dmg = Sg .* rhog - SgOld .* rhogOld;

            % Compute the relative permeability
            [klr, kgr] = constModel.relativePermeabilities(mean(Sl));

            % Derivative of the relative permeability wrt to the saturation
            [dklrdSlm, dkgrdSlm] = constModel.derivativeRelPerm(mean(Sl));

            % Derivative of the relative permeability wrt the pressure
            dklrdPl = -dklrdSlm * dSlmSli * dSldpc;
            dklrdPg =  dklrdSlm * dSlmSli * dSldpc;
            dkgrdPl = -dkgrdSlm * dSlmSli * dSldpc;
            dkgrdPg =  dkgrdSlm * dSlmSli * dSldpc;

            % Initialize the volume of the element
            vol = 0.0;

            % Advective terms
            H = zeros(this.nglp, this.nglp);
            fgrav = zeros(this.nglp, 1);
            for i = 1:this.nIntPoints

                % Shape function matrix
                Np = this.shape.shapeFncMtrx(this.intPoint(i).X);
               
                % Compute the B matrix at the int. point and the detJ
                [Bp, detJ] = this.shape.dNdxMatrix(this.node,this.intPoint(i).X);

                % Get porous media and fluid parameters
                K = this.mat.porousMedia.intrinsicPermeabilityMatrix();

                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;
                if this.isAxisSymmetric
                    c = c * this.shape.axisSymmetricFactor(Np,this.node);
                end

                % Compute permeability matrix
                H = H + Bp' * K * Bp * c;
                
                % Gravity force
                fgrav = fgrav + Bp' * K * grav * c;

                % Compute the element volume
                vol = vol + c;
            end

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

            % Compute the gravity forces
            fel = zeros(this.nnd_el,1);
            feg = zeros(this.nnd_el,1);
            if (this.gravityOn)
                fel = (klr / mul) * mean(rhol) * fgrav;
                feg = (kgr / mug) * mean(rhog) * fgrav;
                Hll = Hll - (1.0/mul) * fgrav * (mean(rhol) * dklrdPl + (klr * dSlmSli / Klb) * rhol)';
                Hlg = Hlg - (1.0/mul) * fgrav * (mean(rhol) * dklrdPg)';
                Hgl = Hgl - (1.0/mug) * fgrav * (mean(rhog) * dkgrdPl)';
                Hgg = Hgg - (1.0/mug) * fgrav * (mean(rhog) * dkgrdPg + (kgr * dSlmSli / Kgb) * rhog)';
            end

            % Storage terms
            masscoeff = phi * (vol / this.nnd_el) / this.DTime;
            fil = fil + (Dml ./ rhol) * masscoeff;
            fig = fig + (Dmg ./ rhog) * masscoeff;
            Cll = ((SlOld .* rholOld ./ rhol) / Klb - dSldpc) * masscoeff;
            Clg = (dSldpc) * masscoeff;
            Cgg = ((SgOld .* rhogOld ./ rhog) / Kgb - dSldpc) * masscoeff;

            % Jacobian matrix
            dfidu = [Hll , Hlg; Hgl , Hgg];

            % Add terms associated with the mass storage
            dfidu = dfidu + [diag(Cll), diag(Clg); diag(Clg), diag(Cgg) ];

            % Assemble element internal force vector
            fi = [fil; fig];

            % Assemble element external force vector
            fe = [fel; feg];

            % Initialize matrices that are not being used
            Ce = zeros(2*this.nglp, 2*this.nglp);
            Ke = zeros(2*this.nglp, 2*this.nglp);
            
        end

        % -----------------------------------------------------------------
        % Compute the permeability tensors
        function [kll, klg, kgl, kgg] = permeabilityTensors(~,ip,pg,pc,Sl)
             [kll, klg, kgl, kgg] = ip.constitutiveMdl.permeabilityMtrcs(Sl,pg-pc,pg);
        end

        % -----------------------------------------------------------------
        % Compute the compressibility coefficients
        function [cll, clg, cgl, cgg] = compressibilityCoeffs(~,ip,pg,pc,Sl)
             [cll, clg, cgl, cgg] =  ip.constitutiveMdl.compressibilityCoeffs(Sl,pg-pc,pg);
        end

        %------------------------------------------------------------------
        % Compute the lumped mass matrices
        function [Sll,Slg,Sgl,Sgg] = lumpedCompressibilityMatrices(this, pc, pg, vol)

            % Shape function matrix
            Np = this.shape.shapeFncMtrx([0.0,0.0]);

            % Pressure values at the integration point
            pcIP = Np * pc;
            pgIP = Np * pg;

            % Compute the saturation degree at the integration point
            Sl = this.intPoint(1).constitutiveMdl.saturationDegree(pcIP);

            % Get compressibility coefficients
            [cll, clg, cgl, cgg] = this.compressibilityCoeffs(this.intPoint(1),pgIP,pcIP,Sl);

            % Mass distribution factor
            factor = vol / this.nnd_el;

            % Compressibility matrices
            Sll = cll * factor * eye(this.nglp,this.nglp);
            Slg = clg * factor * eye(this.nglp,this.nglp);
            Sgl = cgl * factor * eye(this.nglp,this.nglp);
            Sgg = cgg * factor * eye(this.nglp,this.nglp);

        end

        %------------------------------------------------------------------
        % Add contribution of the gravity forces to the external force vector
        function [fel,feg] = addGravityForces(this,fel,feg,Bp,kl,kg,pl,pg,c)

            % Get gravity vector
            grav = this.g * this.mat.porousMedia.b;

            % Get fluid densities
            rhol = this.mat.liquidFluid.getDensity();
            rhog = this.mat.gasFluid.getDensity();

            % Compute the contribution of the gravitational forces
            fel = fel + Bp' * kl * rhol * grav * c;
            feg = feg + Bp' * kg * rhog * grav * c;
            
        end

        %------------------------------------------------------------------
        % Function to get the nodal values of the liquid pressure
        function pl = getNodalLiquidPressure(this)
            pl = this.ue(1:this.nglp);
        end

        %------------------------------------------------------------------
        % Function to get the old nodal values of the liquid pressure
        function plOld = getOldNodalLiquidPressure(this)
            plOld = this.ueOld(1:this.nglp);
        end

        %------------------------------------------------------------------
        % Function to get the nodal values of the gas pressure
        function pg = getNodalGasPressure(this)
            pg = this.ue(1+this.nglp:end);
        end

        %------------------------------------------------------------------
        % Function to get the nodal values of the gas pressure
        function pgOld = getOldNodalGasPressure(this)
            pgOld = this.ueOld(1+this.nglp:end);
        end

        %------------------------------------------------------------------
        % Function to get the nodal values of the capillary pressure
        function pc = getNodalCapillaryPressure(this)
            pl = this.getNodalLiquidPressure();
            pg = this.getNodalGasPressure();
            pc = pg - pl;
        end

        %------------------------------------------------------------------
        % Function to get the nodal values of the capillary pressure
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
