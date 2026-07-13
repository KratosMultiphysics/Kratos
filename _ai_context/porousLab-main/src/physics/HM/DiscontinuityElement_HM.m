%% DiscontinuityElement_HM Class
%
%% Methods
% 
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class Definition
classdef DiscontinuityElement_HM < DiscontinuityElementConductive_HM    
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        ndof_jump = 1;      % Pressure jump dofs
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = DiscontinuityElement_HM(node, mat)
            this = this@DiscontinuityElementConductive_HM(node, mat)
            this.ndof = this.ndof_u + this.ndof_jump + this.ndof_int;
        end
    end

    %% Public methods
    methods

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
        function [fidi, Kddi, Qadi, Hddi, Sddi, Lcci, Lcji, Lcdi, Ljci, Ljji, Ljdi, Ldci, Ldji, Lddi, fi, fe, dfidu] = elementData(this, ae, celem, id)
            
            % Declare output matrices that won't be used
            fi = []; fe = []; dfidu = [];

            % Initialize the matrices for the numerical integration
            fidi = zeros(this.ndof_u, 1);
            Kddi = zeros(this.ndof_u, this.ndof_u);
            Qadi = zeros(this.ndof_u, this.ndof_int);
            Hddi = zeros(this.ndof_int,this.ndof_int);
            Sddi = zeros(this.ndof_int,this.ndof_int);
            Lcci = zeros(celem.nglp,celem.nglp);
            Lcji = zeros(celem.nglp,this.ndof_jump);
            Lcdi = zeros(celem.nglp,this.ndof_int);
            Ljji = zeros(this.ndof_jump,this.ndof_jump);
            Ljdi = zeros(this.ndof_jump,this.ndof_int);
            Lddi = zeros(this.ndof_int,this.ndof_int);

            % Get the discontinuity reference point
            Xr = this.referencePoint();

            % Get the discontinuity tangential vector
            m = this.tangentialVector();

            % Get the normal vector
            nd = this.normalVector();

            % Initialize output matrices
            for i = 1:this.nIntPoints

                % Get the shape function matrix
                N  = this.shape.shapeFnc(this.intPoint(i).X);

                % Cartesian coordinates of the integration point
                Xcar = this.shape.coordNaturalToCartesian(this.node,this.intPoint(i).X);

                % Get the shape function matrix
                Na = this.displacementJumpInterpolation(Xcar,Xr,m);

                % Compute the strain vector
                this.intPoint(i).strain = Na * ae;

                % Compute the stress vector and the constitutive matrix
                [td,Td] = this.intPoint(i).mechanicalLaw();

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
                kh = this.intPoint(i).constitutiveMdl.longitudinalPermeability(this.intPoint(i));

                % Get compressibility coefficient
                comp = this.intPoint(i).constitutiveMdl.compressibility();

                % Get leak-offs
                lt = this.intPoint(i).constitutiveMdl.leakoff;
                lb = this.intPoint(i).constitutiveMdl.leakoff;

                % Numerical integration term. The determinant is ld/2.
                c = detJ * this.intPoint(i).w * this.t;

                % Compute the stiffness sub-matrix
                Kddi = Kddi + Na' * Td * Na * c;

                % Compute the internal force vector
                fidi = fidi + Na' * td * c;

                % Compute the hydromechanical coupling
                Qadi = Qadi + Na' * [0;1] * N * c;

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