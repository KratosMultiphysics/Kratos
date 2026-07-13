%% DiscontinuityElement_H2 Class
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
classdef DiscontinuityElement_H2 < DiscontinuityElement    
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        ndof_int  = 2;      % Discontinuity internal pressure dofs of each phase
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = DiscontinuityElement_H2(node, mat)
            this = this@DiscontinuityElement(node, mat)
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
                constModel = MaterialDiscontinuity_H2(this.mat);
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
        function [Ke, Ce, fi, fe, dfidu] = elementData(this, pl, pg, plOld, pgOld, g, gravityOn, DTime)

            % Get constitutive model
            constModel = this.intPoint(1).constitutiveMdl;

            % Discontinuity geometric parameters
            md = this.tangentialVector();
            ld = this.ld();

            % Get gravity acceleration at the tangential direction
            grav = g * md' * [0.0; -1.0];

            % Get the fluid viscosity
            mul = this.mat.liquidFluid.mu;
            mug = this.mat.gasFluid.mu;

            % Get the fluids bulk modulus
            Klb  = constModel.liquidFluid.K;
            Kgb  = constModel.gasFluid.K;
            
            % Vector of the nodal pore-pressure dofs
            pc = pg - pl;

            % Vector with the old nodal dofs
            pcOld = pgOld - plOld;

            % Fill nodal state variables
            Sl      = zeros(2,1);
            SlOld   = zeros(2,1);
            dSldpc  = zeros(2,1);
            rhol    = zeros(2,1);
            rhog    = zeros(2,1);
            rholOld = zeros(2,1);
            rhogOld = zeros(2,1);
            for i = 1:2
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
            dSlmSli = 0.5;

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

            % Get the longitudina permeability
            K = constModel.longitudinalPermeability();

            % Get the aperture
            aperture = constModel.initialAperture;

            % Get the discontinuity porosity
            phid = constModel.porosity;

            % Advective terms
            H = zeros(this.ndof_int, this.ndof_int);
            fgrav = zeros(this.ndof_int, 1);
            for i = 1:this.nIntPoints
               
                % Compute the B matrix at the int. point and the detJ
                [Bp, detJ] = this.shape.dNdxMatrix(this.node,this.intPoint(i).X);

                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;

                % Compute permeability matrix
                H = H + Bp' * K * aperture * Bp * c;
                
                % Gravity force
                fgrav = fgrav + Bp' * K * aperture * grav * c;  % check
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
            fel = zeros(this.ndof_int,1);
            feg = zeros(this.ndof_int,1);
            if (gravityOn)
                fel = (klr / mul) * mean(rhol) * fgrav;
                feg = (kgr / mug) * mean(rhog) * fgrav;
                Hll = Hll - (1.0/mul) * fgrav * (mean(rhol) * dklrdPl + (klr * dSlmSli / Klb) * rhol)';
                Hlg = Hlg - (1.0/mul) * fgrav * (mean(rhol) * dklrdPg)';
                Hgl = Hgl - (1.0/mug) * fgrav * (mean(rhog) * dkgrdPl)';
                Hgg = Hgg - (1.0/mug) * fgrav * (mean(rhog) * dkgrdPg + (kgr * dSlmSli / Kgb) * rhog)';
            end

            % Storage terms
            masscoeff = phid * aperture * (ld / 2.0) / DTime;
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
            Ce = zeros(2*this.ndof_int, 2*this.ndof_int);
            Ke = zeros(2*this.ndof_int, 2*this.ndof_int);
        end
    end
end