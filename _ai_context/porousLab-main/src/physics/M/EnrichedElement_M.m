%% EnrichedElement_M Class
% This class defines an enriched mechanical finite element that extends 
% the functionality of the _RegularElement_M_ class. It incorporates 
% additional degrees of freedom to handle discontinuities and enrichment 
% modes such as stretching and relative rotation. 
%
%% Methods
% * *elementData*: Computes the element data (stiffness matrix, damping 
%                  matrix, internal force vector, external force vector, 
%                  and derivative of internal force vector) based on 
%                  whether the element has discontinuities.
% * *enrichedElementData*: Computes the enriched element data using static 
%                          condensation of the enrichment degrees of 
%                          freedom.
% * *solveLocalEq*: Solves the local equilibrium equation iteratively 
%                   using the Newton-Raphson method.
% * *fillElementSubData*: Fills the sub-matrices and vectors for the 
%                         enriched element based on the current enriched 
%                         degrees of freedom.
% * *getDiscontinuitiesData*: Computes the stiffness matrix and force 
%                             vector contributions from the 
%                             discontinuities.
% * *getNumberEnrichedDofs*: Returns the total number of enriched degrees 
%                            of freedom.
% * *getNumberOfDiscontinuities*: Returns the number of discontinuities 
%                                 associated with the element.
% * *getNumberOfDofPerDiscontinuity*: Returns the number of degrees of 
%                                     freedom per discontinuity, 
%                                     considering the enabled enrichment 
%                                     modes.
% * *addDiscontinuitySegment*: Adds a discontinuity segment to the element.
% * *kinematicEnrichment*: Computes the kinematic enrichment matrix for 
%                          the element based on the discontinuities and 
%                          enrichment modes.
% 
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00: Initial version (January 2024).
%
%% Class definition
classdef EnrichedElement_M < RegularElement_M    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        discontinuity                = [];
        addTangentialStretchingMode  = false;
        addNormalStretchingMode      = false;
        addRelRotationMode           = false;
        condenseEnrDofs              = true;
        symmetricForm                = true;
        stressIntCoeff               = [];
        useNodalEnrDofs              = false;
        nNodalEnrDofs                = 4;
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = EnrichedElement_M(node, elem, t, ...
                mat, intOrder, glu, massLumping, lumpStrategy, ...
                isAxisSymmetric,isPlaneStress, ...
                addRelRotationMode,addTangentialStretchingMode, ...
                addNormalStretchingMode, condenseEnrDofs,...
                subDivInt, symmetricForm, useNodalEnrDofs)
            this = this@RegularElement_M(node, elem, t, ...
                mat, intOrder, glu, massLumping, lumpStrategy, ...
                isAxisSymmetric,isPlaneStress);
            this.addTangentialStretchingMode  = addTangentialStretchingMode;
            this.addNormalStretchingMode      = addNormalStretchingMode;
            this.addRelRotationMode           = addRelRotationMode;
            this.condenseEnrDofs              = condenseEnrDofs;
            this.subDivInt                    = subDivInt;
            this.symmetricForm                = symmetricForm;
            this.useNodalEnrDofs              = useNodalEnrDofs;
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Computes the element data for the current element based on whether
        % the element contains a discontinuity or not.
        % 
        % Outputs:
        %   Ke    - Element stiffness matrix.
        %   Ce    - Element damping matrix.
        %   fi    - Internal force vector.
        %   fe    - External force vector.
        %   dfidu - Derivative of internal force with respect to 
        %           displacement.
        function [Ke, Ce, fi, fe, dfidu] = elementData(this)

           if isempty(this.discontinuity)
               [Ke, Ce, fi, fe, dfidu] = elementData@RegularElement_M(this);
           else
               [Ke, Ce, fi, fe, dfidu] = enrichedElementData(this);
           end
            
        end

        % -----------------------------------------------------------------
        % Update state variables.
        function updateStateVar(this)
            updateStateVar@RegularElement_M(this);
            % Loop through the discontinuities
            nDiscontinuities = this.getNumberOfDiscontinuities();
            if nDiscontinuities == 0, return, end
            for i = 1:nDiscontinuities
                this.discontinuity(i).updateStateVar();
            end
        end

        %------------------------------------------------------------------
        % Function to reset the displacements and strains
        function udofs = resetDisplacements(this)

            udofs = this.gle;

            resetDisplacements@RegularElement_M(this);

            % Loop through the discontinuities
            nDiscontinuities = this.getNumberOfDiscontinuities();
            if nDiscontinuities == 0, return, end
            for i = 1:nDiscontinuities
                this.discontinuity(i).resetDisplacements();
            end
        end

        %------------------------------------------------------------------
        % Computes the enriched element data for the finite element.
        % 
        % Outputs:
        %   Ke    - Element stiffness matrix.
        %   Ce    - Element damping matrix.
        %   fi    - Internal force vector.
        %   fe    - External force vector.
        %   dfidu - Derivative of internal force with respect to 
        %           displacement.
        function [Ke, Ce, fi, fe, dfidu] = enrichedElementData(this)

            % Initialize the matrices and vectors that will be returned
            Ce = zeros(this.ngle, this.ngle);
            Ke = zeros(this.ngle, this.ngle);

            if this.condenseEnrDofs
                % Compute the sub-matrices
                [Kuu, Kua, Kau, Kaa, fiu, fia, feu] = this.solveLocalEq();
    
                % Static condensation of the enrichment dofs
                dfidu = Kuu - Kua * (Kaa\Kau);
                fi    = fiu - Kua * (Kaa\fia);
                fe    = feu;
            else
                % Get enrichment dofs
                ae = this.ue(this.nglu+1:end);

                % Compute sub-matrices
                [Kuu, Kua, Kau, Kaa, fiu, fia, feu, fea] = this.fillElementSubData(ae);

                % Assemble element matrices
                dfidu = [ Kuu , Kua;
                          Kau , Kaa];
                fi    = [ fiu; fia];
                fe    = [ feu; fea];

            end

        end

        %------------------------------------------------------------------
        %  Solves the local equilibrium equation for an enriched element
        function [Kuu, Kua, Kau, Kaa, fiu, fia, feu, fea] = solveLocalEq(this)

            % Initialize the enrichment dofs vector
            nEnrDofs = this.getNumberEnrichedDofs();
            ae = zeros(nEnrDofs,1);

            % Define Newton-Raphson parameter
            conv    = false;
            tol     = 1.0e-5;
            maxIter = 10;

            % Iterative solution of the local equilibrium equation
            for i = 1:maxIter

                % Compute sub-matrices
                [Kuu, Kua, Kau, Kaa, fiu, fia, feu, fea] = this.fillElementSubData(ae);

                % Check convergence
                if (norm(fia) < tol)
                    conv = true;
                    break
                end

                % Solve iterative equation
                dae = -Kaa\fia;

                % Update enriched dofs
                ae = ae + dae;
                    
            end

            % Stop analysis process if solution did not converge
            if (conv == false)
                disp('LOCAL EQUILIBRIUM EQUATION FAIL TO CONVERGE');
                error('Local equilibrium equation did not converge');
            end

        end

        %------------------------------------------------------------------
        % Computes and assembles the sub-matrices and sub-vectors for an
        % enriched finite element
        function [Kuu, Kua, Kau, Kaa, fiu, fia, feu, fea] = fillElementSubData(this,ae)

            % Get the number of dofs of each type
            nEnrDofs = this.getNumberEnrichedDofs();
            nRegDofs = this.nglu;

            % Initialize the sub-matrices
            Kuu = zeros(nRegDofs,nRegDofs);
            Kua = zeros(nRegDofs,nEnrDofs);
            Kau = zeros(nEnrDofs,nRegDofs);
            Kaa = zeros(nEnrDofs,nEnrDofs);

            % Initialize the sub-vectors
            fiu = zeros(nRegDofs,1);
            fia = zeros(nEnrDofs,1);

            % Initialize external force vector
            feu = zeros(nRegDofs, 1);
            fea = zeros(nEnrDofs, 1);

            % Vector of the nodal dofs
            u  = this.getNodalDisplacement();

            % Numerical integration of the sub-matrices
            for i = 1:this.nIntPoints
               
                % Compute the B matrix at the int. point and the detJ
                [dNdx, detJ] = this.shape.dNdxMatrix(this.node,this.intPoint(i).X);

                % Assemble the B-matrix for the mechanical part
                Bu = this.BMatrix(dNdx);

                % Get kinematic enriched matrix
                Gr = this.kinematicEnrichment(Bu,this.intPoint(i).X);

                % Get the static enriched matrix
                if this.symmetricForm
                    Gv = Gr;
                else
                    Gv = this.equilibriumEnrichment(this.intPoint(i).X);
                end

                % Compute the strain vector
                this.intPoint(i).strain = Bu * u + Gr * ae;

                % Compute the stress vector and the constitutive matrix
                [stress,Duu] = this.intPoint(i).mechanicalLaw();
        
                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;
                if this.isAxisSymmetric
                    c = c * this.shape.axisSymmetricFactor(Np,this.node);
                end
                
                % Compute the stiffness sub-matrix
                Kuu = Kuu + Bu' * Duu * Bu * c;
                Kua = Kua + Bu' * Duu * Gr * c;
                Kau = Kau + Gv' * Duu * Bu * c;
                Kaa = Kaa + Gv' * Duu * Gr * c;

                % Internal force vector
                fiu = fiu + Bu' * stress * c;
                fia = fia + Gv' * stress * c;
                
                % Compute the gravity forces
                if (this.gravityOn)
                    feu = this.addGravityForces(feu,this.intPoint(i).X,c);
                end
            end

            % Add the contribution from the discontinuities
            [Kd,fd] = this.getDiscontinuitiesData(ae);
            fia = fia + fd;
            Kaa = Kaa + Kd;
        end

        %------------------------------------------------------------------
        % Computes the stiffness matrix and force vector contributions
        % from discontinuities in the enriched element
        function [Kd,fd] = getDiscontinuitiesData(this,ae)

            nEnrDofs          = this.getNumberEnrichedDofs();
            nDiscontinuities  = this.getNumberOfDiscontinuities();
            nLocalDofDiscontinuity = this.getNumberOfDofPerDiscontinuity();
            nDofDiscontinuity = nLocalDofDiscontinuity;
            if this.useNodalEnrDofs
                nDofDiscontinuity = this.nNodalEnrDofs;
            end

            % Initialize the output data 
            Kd = zeros(nEnrDofs,nEnrDofs);
            fd = zeros(nEnrDofs,1);

            % Loop through the discontinuities
            for i = 1:nDiscontinuities

                % Dofs associated with this discontinuity segment
                dofs = nDofDiscontinuity*(i-1)+1 : nDofDiscontinuity*i;

                % Get the discontinuity data
                [Kdi,~,fdi,~,~] = this.discontinuity(i).elementData(ae(dofs));

                % Assemble the contribution of this discontinuity
                Kd(dofs,dofs) = Kdi;
                fd(dofs)      = fdi;
                
            end
        end

        %------------------------------------------------------------------
        % Compute the forces due to the pore-pressure field
        function fe = porePressureForce(this, pe)

            if isempty(this.discontinuity)
               fe = porePressureForce@RegularElement_M(this, pe);
               return;
            end

            % Get the number of dofs of each type
            nEnrDofs = this.getNumberEnrichedDofs();
            nRegDofs = this.nglu;

            % Initialize hydro-mechanical coupling forces
            feu = zeros(nRegDofs,1);
            fea = zeros(nEnrDofs,1);

            % Identity vector
            m = [1;1;1;0];

            % Numerical integration of Q
            for i = 1:this.nIntPoints

                % Shape function vector
                N = this.shape.shapeFncMtrx(this.intPoint(i).X);

                % Compute the B matrix at the int. point and the detJ
                [dNdx, detJ] = this.shape.dNdxMatrix(this.node,this.intPoint(i).X);

                % Assemble the B-matrix for the mechanical part
                Bu = this.BMatrix(dNdx,N);

                % Get the static enriched matrix
                if this.symmetricForm
                    Gv = this.kinematicEnrichment(Bu,this.intPoint(i).X);
                else
                    Gv = this.equilibriumEnrichment(this.intPoint(i).X);
                end

                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;
                if this.isAxisSymmetric
                    c = c * this.shape.axisSymmetricFactor(N,this.node);
                end

                % Evaluate pore-pressure at the integration point
                p = this.porePressureField(N, this.intPoint(i).X, pe);

                % Compute the forces
                feu = feu + Bu' * m * p * c;
                fea = fea + Gv' * m * p * c;

            end

            % Loop through the discontinuities
            nDiscontinuities  = this.getNumberOfDiscontinuities();
            nDofDiscontinuity = this.getNumberOfDofPerDiscontinuity();
            for i = 1:nDiscontinuities

                % Initialize the discontinuity pore-pressure vector
                pd = [0; 0];
                for j = 1:2
                    % Cartesian coordinate of node i of the discontinuity
                    X = this.discontinuity(i).node(j,:);
                    % Natural coordinates of this node
                    Xn = this.shape.coordCartesianToNatural(this.node,X);
                    % Shape function
                    N = this.shape.shapeFncMtrx(Xn);
                    % Function phi
                    phi = this.auxiliaryfncPhi(i, N);
                    % Pore-pressure value at each side of the discontinuity
                    ptop = N * pe + (1.0-phi)*this.discontinuity(i).DP;
                    pbot = N * pe - phi*this.discontinuity(i).DP;
                    % Pore-pressure at the discontinuity
                    pd(j) = mean([ptop, pbot]);
                end

                % Get the discontinuity data
                feai = this.discontinuity(i).porePressureForce(pd);

                % Dofs associated with this discontinuity segment
                dofs = nDofDiscontinuity*(i-1)+1 : nDofDiscontinuity*i;

                % Assemble the contribution of this discontinuity
                fea(dofs,1) = fea(dofs,1) + feai;
                
            end

            % Assemble vector
            fe = [feu; fea];

        end

        %------------------------------------------------------------------
        % Evaluates the pore-pressure field
        function p = porePressureField(this, N, Xn, pe)

            % Continuous part
            p = N * pe;

            % Cartesian coordinate of the integration point
            X = this.shape.coordNaturalToCartesian(this.node, Xn);

            % Discontinuous part
            nDiscontinuities = this.getNumberOfDiscontinuities();
            for i = 1:nDiscontinuities
                % Function phi
                phi = this.auxiliaryfncPhi(i, N);
                % Heaviside function
                h = this.discontinuity(i).heaviside(X);
                % Add discontinuous part from discontinuity i
                p = p + (h-phi) * this.discontinuity(i).DP;
            end
        end

        %------------------------------------------------------------------
        % Gets the number of enriched degrees of freedom
        function nEnrDof = getNumberEnrichedDofs(this)
            nEnrDof = this.getNumberOfDiscontinuities();
            nLocalDofDiscontinuity = this.getNumberOfDofPerDiscontinuity();
            nDofDiscontinuity = nLocalDofDiscontinuity;
            if this.useNodalEnrDofs
                nDofDiscontinuity = this.nNodalEnrDofs;
            end
            nEnrDof = nEnrDof * nDofDiscontinuity;
        end

        %------------------------------------------------------------------
        % Obtain the number of discontinuities
        function n = getNumberOfDiscontinuities(this)
            n = size(this.discontinuity,1);
        end

        %------------------------------------------------------------------
        % Calculates the number of degrees of freedom per discontinuity for
        % the enriched element
        function n = getNumberOfDofPerDiscontinuity(this)
            n = 2;  
            if this.addTangentialStretchingMode
                n = n + 1;
            end
            if this.addNormalStretchingMode
                n = n + 1;
            end
            if this.addRelRotationMode
                n = n + 1;
            end
        end

        %------------------------------------------------------------------
        % Displacement jump order
        function n = displacementJumpOrder(this)
            n = 0;
            if (this.addTangentialStretchingMode || this.addNormalStretchingMode || this.addRelRotationMode) 
                n = 1;
            end
        end

        %------------------------------------------------------------------
        % Adds a discontinuity segment to the element
        function addDiscontinuitySegment(this,dseg)
            this.discontinuity = [this.discontinuity; dseg];
        end

        %------------------------------------------------------------------
        % Adds the discontinuities dofs to the element dof vector
        function addEnrichmentToDofVector(this)
            nDiscontinuities = this.getNumberOfDiscontinuities();
            for i = 1:nDiscontinuities
                this.gle = [this.gle, this.discontinuity(i).dof];
            end
            this.ngle = length(this.gle);
        end

        %------------------------------------------------------------------
        % Computed the kinematic enrichment matrix for an enriched finite
        % element.
        % It calculated the enrichment matrix by considering the
        % contributions of discontinuities in the element. The enrichment
        % includes translation, stretching and relative rotation modes
        % depending on the configuration
        function Gc = kinematicEnrichment(this, Bu,Xn) 
            nDiscontinuities  = this.getNumberOfDiscontinuities();
            nLocalDofDiscontinuity = this.getNumberOfDofPerDiscontinuity();
            nDofDiscontinuity = nLocalDofDiscontinuity;
            if this.useNodalEnrDofs
                nDofDiscontinuity = this.nNodalEnrDofs;
            end
            Gc = zeros(4,nDofDiscontinuity * nDiscontinuities);
            for i = 1:nDiscontinuities    
                Gci = zeros(4,nLocalDofDiscontinuity);
                % Get the discontinuity orientation vectors
                m = this.discontinuity(i).tangentialVector();
                n = this.discontinuity(i).normalVector();
                % Get the discontinuity reference point
                Xr = this.discontinuity(i).referencePoint();
                for j = 1:this.nnd_el
                    Xj = this.node(j,:);
                    h = this.discontinuity(i).heaviside(Xj);
                    if (h > 0.0)
                        % Columns of the B-matrix associated with this node
                        Buj = Bu(:, 2*(j-1) + 1 : 2*j);
                        % Add translation modes
                        Gci(:,1) = Gci(:,1) - Buj * m;
                        Gci(:,2) = Gci(:,2) - Buj * n;
                        % Add stretching mode
                        c = 3;
                        if this.addTangentialStretchingMode
                            Gci(:,c) = Gci(:,c) - Buj * (m * m') * (Xj' - Xr');
                            c = c + 1;
                        end
                        if this.addNormalStretchingMode
                            Gci(:,c) = Gci(:,c) - Buj * (n * n') * (Xj' - Xr');
                            c = c + 1;
                        end
                        % Add relative rotation mode
                        if this.addRelRotationMode
                            mn = (n * m') - (m * n');
                            Gci(:,c) = Gci(:,c) - Buj * mn * (Xj' - Xr');
                        end
                    end
                end
                % Add stretching mode
                if (this.addTangentialStretchingMode || this.addNormalStretchingMode)
                    % Cartesian coordinate of the int. point
                    X = this.shape.coordNaturalToCartesian(this.node,Xn);
                    % Heaviside function
                    h = this.discontinuity(i).heaviside(X);
                    % Fill matrix Gci
                    c = 3;
                    if this.addTangentialStretchingMode
                        M = [m(1),0;0,m(2);0,0;m(2),m(1)];
                        Gci(:,c) = Gci(:,c) + h * M * m;
                        c = c + 1;
                    end
                    if this.addNormalStretchingMode
                        N = [n(1),0;0,n(2);0,0;n(2),n(1)];
                        Gci(:,c) = Gci(:,c) + h * N * n;
                    end                    
                end
                % Assemble the matrix associated with discontinuity i
                cols = nDofDiscontinuity*(i-1)+1 : nDofDiscontinuity*i;
                % Get transformation matrix from local to global
                T = this.discontinuity(i).getDofTransformationMtrx();
                Gc(:,cols) = Gci * T;
            end
        end

        %------------------------------------------------------------------
        % Compute the enrichment auxiliary function
        function phi = auxiliaryfncPhi(this, id, N) 
            phi = 0.0;
            for j = 1:this.nnd_el
                Xj = this.node(j,:);
                hj = this.discontinuity(id).heaviside(Xj);
                if (hj > 0.0)
                    phi = phi + N(:, j);
                end
            end
        end

        %------------------------------------------------------------------
        % Compute the enrichment auxiliary function
        function Gv = equilibriumEnrichment(this, Xn)
            
            % Initialize variables
            nDiscontinuities  = this.getNumberOfDiscontinuities();
            nLocalDofDiscontinuity = this.getNumberOfDofPerDiscontinuity();
            nDofDiscontinuity = nLocalDofDiscontinuity;
            if this.useNodalEnrDofs
                nDofDiscontinuity = this.nNodalEnrDofs;
            end
            jumpOrder = this.displacementJumpOrder();
            Gv = zeros(4,nDofDiscontinuity * nDiscontinuities);

            % Get the stress interpolation coefficients
            c = this.getStressInterpCoefficients();

            % Cartesian coordinate of the int. point
            X = this.shape.coordNaturalToCartesian(this.node,Xn);

            % Stress interpolation polynomial
            p = this.shape.polynomialStress(X);

            % Evaluate the polynomial
            g = c'*p;

            % Loop through the discontinuities
            for i = 1:nDiscontinuities
                
                % Get discontinuity data
                m = this.discontinuity(i).tangentialVector();
                n = this.discontinuity(i).normalVector();
                P = this.discontinuity(i).projectionMatrix();
                
                % Get matrix associated with discontinuity i
                if (jumpOrder == 0)
                    Gi = - g(1,i) * P * [ m , n];
                elseif (jumpOrder == 1)
                    if ((this.addTangentialStretchingMode == true) && (this.addRelRotationMode == false))
                        Gi = - P * [ g(1,i)*m , g(1,i)*n, g(2,i)*m];
                    elseif ((this.addTangentialStretchingMode == false) && (this.addRelRotationMode == true))
                        Gi = - P * [ g(1,i)*m , g(1,i)*n, g(2,i)*n];
                    else
                        Gi = - P * [ g(1,i)*m , g(1,i)*n, g(2,i)*m, g(2,i)*n];
                    end
                end

                % Assemble the matrix associated with discontinuity i
                cols = nDofDiscontinuity*(i-1)+1 : nDofDiscontinuity*i;

                % Get transformation matrix from local to global
                T = this.discontinuity(i).getDofTransformationMtrx();

                Gv(:,cols) = Gi * T;
            end

        end

        %------------------------------------------------------------------
        % Get the stress interpolation coefficients
        function c = getStressInterpCoefficients(this)
            if isempty(this.stressIntCoeff)
                this.stressIntCoeff = this.computeStressIntCoeffs();
            end
            c = this.stressIntCoeff;
        end

        %------------------------------------------------------------------
        % Compute the stress interpolation coefficients
        function c = computeStressIntCoeffs(this)
            
            % Gramm matrix
            H = this.grammMatrix();

            % Integral of the stresses along the discontinuities
            S = this.dSetIntegralPolynomialStress();

            % Compute coefficients
            c = H\S;
        end

        %------------------------------------------------------------------
        % Compute the stress interpolation coefficients
        function H = grammMatrix(this)

            % Initialize variables
            n = this.shape.dimPolynomialStressInterp();
            H = zeros(n,n);

            for i = 1:this.nIntPoints
                % Numerical int. coefficient
                detJ = this.shape.detJacobian(this.node,this.intPoint(i).X);
                c = this.intPoint(i).w * detJ * this.t;

                % Cartesian coordinate of the int. point
                X = this.shape.coordNaturalToCartesian(this.node,this.intPoint(i).X);

                % Stress interpolation polynomial
                p = this.shape.polynomialStress(X);

                % Gram matrix
                H = H + (p * p') * c;
            end
        end

        %------------------------------------------------------------------
        % Compute the stress interpolation coefficients
        function S = dSetIntegralPolynomialStress(this)

            % Initialize variables
            dimPolyStress    = this.shape.dimPolynomialStressInterp();
            nDiscontinuities = this.getNumberOfDiscontinuities();
            jumpOrder        = this.displacementJumpOrder();

            % Initialize matrix
            S = zeros(dimPolyStress,nDiscontinuities*(jumpOrder + 1));

            % Loop through the discontinuities
            for i = 1:nDiscontinuities

                % Evaluate integral of discontinuity i
                Si = this.discontinuity(i).intPolynomialStressIntp(this);

                % Assemble
                cols = ((i-1) * (jumpOrder + 1) + 1):(i * (jumpOrder + 1));
                S(:,cols) = Si;
            end
        end

    end
end
