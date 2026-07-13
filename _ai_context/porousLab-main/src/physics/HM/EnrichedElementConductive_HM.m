%% EnrichedElementConductive_HM Class
% This class extends the _RegularElement_HM_ class to define a conductive
% hydro-mechanical finite element with displacement enrichments for
% discontinuities. It provides methods to compute enriched element data,
% manage discontinuity segments, and calculate enriched degrees of freedom.
%
%% Methods
% * *elementData*: Computes the element data (stiffness matrix, damping 
%                  matrix, internal force vector, external force vector, 
%                  and derivative of internal force vector) based on 
%                  whether the element has discontinuities.
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
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class definition
classdef EnrichedElementConductive_HM < RegularElement_HM   
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        discontinuity = [];
        addTangentialStretchingMode  = false;
        addNormalStretchingMode      = false;
        addRelRotationMode           = false;
        symmetricForm                = true;
        stressIntCoeff               = [];
        condenseInternalPressure     = false;
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = EnrichedElementConductive_HM(node, elem, t, ...
                mat, intOrder, glu, glp, massLumping, lumpStrategy, ...
                isAxisSymmetric,isPlaneStress, ...
                addRelRotationMode,addTangentialStretchingMode, ...
                addNormalStretchingMode,...
                subDivInt, symmetricForm)
            this = this@RegularElement_HM(node, elem, t, ...
                mat, intOrder, glu, glp, massLumping, lumpStrategy, ...
                isAxisSymmetric,isPlaneStress);
            this.addTangentialStretchingMode  = addTangentialStretchingMode;
            this.addNormalStretchingMode      = addNormalStretchingMode;
            this.addRelRotationMode           = addRelRotationMode;
            this.subDivInt                    = subDivInt;
            this.symmetricForm                = symmetricForm;
        end
    end
    
    %% Public methods
    methods

        % -----------------------------------------------------------------
        % Update state variables.
        function updateStateVar(this)
            updateStateVar@RegularElement_HM(this);
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

            udofs = this.glu;

            resetDisplacements@RegularElement_HM(this);

            % Loop through the discontinuities
            nDiscontinuities = this.getNumberOfDiscontinuities();
            if nDiscontinuities == 0, return, end
            for i = 1:nDiscontinuities
                this.discontinuity(i).resetDisplacements();
                udofs = [udofs,this.discontinuity(i).dof];
            end
        end

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
               % Get the continuum contribution
               [Ke, Ce, fi, fe, dfidu] = elementData@RegularElement_HM(this);
           else
               % Compute contribution of the continuum
               [Kuu, Kua, Kau, Kaa, Quc, Qac, Qca, Hcc, Scc, fiu, fia, fic, feu, fea, fec] = this.fillElementSubData();
               
               % Compute contribution of the discontinuity
               [fid, Kdd, Qad, Hdd, Sdd, Tdc] = getDiscontinuitiesData(this);

               % Add contribution of the discontinuity stiffness
               fia = fia + fid;
               Kaa = Kaa + Kdd;

               % Condense the internal discontinuity pressure
               Hcc = Hcc + Tdc' * Hdd * Tdc;
               Scc = Scc + Tdc' * Sdd * Tdc;
               Qac = Qac + Qad * Tdc;
               Qca = Qca + Tdc' * Qad';

               % Number of displacement and pressure jump dofs
               ndofa = size(Kaa,1);

               % Create auxiliary zero sub-matrices
               Ouu = zeros(this.nglu, this.nglu);
               Ouc = zeros(this.nglu, this.nglp);
               Oua = zeros(this.nglu, ndofa);

               Ocu = zeros(this.nglp, this.nglu);
               Occ = zeros(this.nglp, this.nglp);
               Oca = zeros(this.nglp, ndofa);

               Oau = zeros(ndofa, this.nglu);
               Oac = zeros(ndofa, this.nglp);
               Oaa = zeros(ndofa, ndofa);
               
               % Assemble the element matrices
               dfidu = [ Kuu,  Ouc,  Kua
                         Ocu,  Occ,  Oca;
                         Kau,  Oac,  Kaa];

               Ke = [ Ouu,  -Quc,  Oua;
                      Ocu,   Hcc,  Oca;
                      Oau,  -Qac,  Oaa];
               
               Ce = [ Ouu,  Ouc,  Oua;
                      Quc', Scc,  Qca;
                      Oau,  Oac,  Oaa];
               
               fi = [fiu; fic; fia];
               
               fe = [feu; fec; fea];
           
           end
        end

        %------------------------------------------------------------------
        % Computes and assembles the sub-matrices and sub-vectors for an
        % enriched finite element
        function [Kuu, Kua, Kau, Kaa, Quc, Qac, Qca, Hcc, Scc, fiu, fia, fic, feu, fea, fec] = fillElementSubData(this)

            % Number of displacement and pressure jump dofs
            nd = this.getNumberOfDiscontinuities();
            ndofa = nd * this.getNumberOfDisplacementDofPerDiscontinuity();

            % Initialize stiffness sub-matrices
            Kuu = zeros(this.nglu, this.nglu);
            Kua = zeros(this.nglu, ndofa);
            Kau = zeros(ndofa, this.nglu);
            Kaa = zeros(ndofa, ndofa);

            % Initialize the hydro-mechanical coupling sub-matrices
            Quc = zeros(this.nglu, this.nglp);
            Qac = zeros(ndofa, this.nglp);
            Qca = zeros(this.nglp, ndofa);
            
            % Initialize fluid-flow sub-matrices
            Hcc = zeros(this.nglp, this.nglp);

            % Initialize compressibity sub-matrices
            Scc = zeros(this.nglp , this.nglp);

            % Initialize external force vector
            feu = zeros(this.nglu, 1);
            fec = zeros(this.nglp, 1);
            fea = zeros(ndofa, 1);

            % Initialize internal force vector
            fiu = zeros(this.nglu, 1);
            fic = zeros(this.nglp, 1);
            fia = zeros(ndofa, 1);
            
            % Vector of the nodal dofs
            u  = this.getNodalDisplacement();
            pc = this.getNodalPressure();
            a  = this.getDisplacementJumpDofs(ndofa);

            % Initialize 2D identity vector
            m = [1.0 ; 1.0 ; 1.0 ; 0.0];

            % Numerical integration of the sub-matrices
            for i = 1:this.nIntPoints

                % Shape function matrix
                Npc = this.shape.shapeFncMtrx(this.intPoint(i).X);

                % Cartesian coordinates of the integration point
                Xcar = this.shape.coordNaturalToCartesian(this.node,this.intPoint(i).X);

                % Enriched shape function matrix
                Npj = this.enrichedShapeFncMtrx(Npc, Xcar);
               
                % Compute the B matrix at the int. point and the detJ
                [Bp, detJ] = this.shape.dNdxMatrix(this.node,this.intPoint(i).X);

                % Assemble the B-matrix for the mechanical part
                Bu = this.BMatrix(Bp);

                % Get kinematic enriched matrix
                Gr = this.kinematicEnrichment(Bu);

                % Get the static enriched matrix
                if this.symmetricForm
                    Gv = Gr;
                else
                    Gv = this.equilibriumEnrichment(this.intPoint(i).X);
                end

                % Pressure values at the integration point
                pIP = Npc * pc;
        
                % Compute the permeability matrix
                kh = this.intPoint(i).constitutiveMdl.permeabilityTensor();

                % Get compressibility coefficient
                comp = this.intPoint(i).constitutiveMdl.compressibilityCoeff();

                % Get Biot's coefficient
                biot = this.intPoint(i).constitutiveMdl.biotCoeff();

                % Compute the strain vector
                this.intPoint(i).strain = Bu * u + Gr * a;

                % Compute the stress vector and the constitutive matrix
                [stress,Duu] = this.intPoint(i).mechanicalLaw();
        
                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;
                if this.isAxisSymmetric
                    c = c * this.shape.axisSymmetricFactor(Np,this.node);
                end

                % Internal force vector
                fiu = fiu + Bu' * stress * c;
                fia = fia + Gv' * stress * c;

                % Compute the stiffness sub-matrices
                Kuu = Kuu + Bu' * Duu * Bu * c;
                Kua = Kua + Bu' * Duu * Gr * c;
                Kau = Kau + Gv' * Duu * Bu * c;
                Kaa = Kaa + Gv' * Duu * Gr * c;

                % Compute the hydro-mechanical coupling sub-matrices
                Quc = Quc + Bu'  * m  * biot * Npc * c;
                Qac = Qac + Gv'  * m  * biot * Npc * c;
                Qca = Qca + Npc' * m' * biot * Gr  * c;
        
                % Compute permeability sub-matrices
                Hcc = Hcc + Bp' * kh * Bp * c;

                % Compute compressibility matrices
                Scc = Scc + Npc' * comp * Npc * c;
                
                % Compute the gravity forces
                if (this.gravityOn)
                    [feu,fep] = this.addGravityForces(feu,fep,Npc,Bp,kh,pIP,c);
                end
            end

        end

        %------------------------------------------------------------------
        % Computes the stiffness matrix and force vector contributions
        % from discontinuities in the enriched element
        function [fid, Kdd, Qad, Hdd, Sdd, Tdc] = getDiscontinuitiesData(this)

            % Number of displacement and pressure jump dofs
            nDiscontinuities = this.getNumberOfDiscontinuities();
            ndofa_discontinuity = this.getNumberOfDisplacementDofPerDiscontinuity();
            ndofa = nDiscontinuities * ndofa_discontinuity;
            ndofd = nDiscontinuities * 2;

            % Initialize the output data
            fid = zeros(ndofa, 1);
            Kdd = zeros(ndofa, ndofa);
            Qad = zeros(ndofa, ndofd);
            Hdd = zeros(ndofd, ndofd);
            Sdd = zeros(ndofd, ndofd);

            % Condensation matrices
            Tdc = zeros(ndofd,this.nglp);
            
            % Get dof vector
            a = this.getDisplacementJumpDofs(ndofa);

            % Loop through the discontinuities
            for i = 1:nDiscontinuities

                % Displacement jump dofs associated with this discontinuity
                % segment (local numbering)
                dofs_a =  ndofa_discontinuity*(i-1)+1 : ndofa_discontinuity*i;

                % Internal pressure dofs associated with this discontinuity
                % segment (local numbering)
                dofs_d = 2*(i-1)+1 : 2*i;

                % Get the discontinuity data
                [fidi, Kddi, Qadi, Hddi, Sddi,~,~,~] = this.discontinuity(i).elementData(a(dofs_a));

                % Assemble the contribution of this discontinuity
                fid(dofs_a) = fid(dofs_a) + fidi;
                Kdd(dofs_a, dofs_a) = Kdd(dofs_a, dofs_a) + Kddi;
                Qad(dofs_a, dofs_d) = Qad(dofs_a, dofs_d) + Qadi;
                Hdd(dofs_d, dofs_d) = Hdd(dofs_d, dofs_d) + Hddi;
                Sdd(dofs_d, dofs_d) = Sdd(dofs_d, dofs_d) + Sddi;

                % Loop through the nodes of the discontinuity to fill the
                % condensation matrix
                for j = 1:2
                    X = this.discontinuity(i).node(j,:);
                    Xn = this.shape.coordCartesianToNatural(this.node,X);
                    Np = this.shape.shapeFncMtrx(Xn);
                    Tdc(dofs_d(j),:) = Np;
                end
                
            end
        end

        %------------------------------------------------------------------
        % Gets the number of enriched degrees of freedom
        function nEnrDof = getNumberEnrichedDofs(this)
            nEnrDof = this.getNumberOfDiscontinuities();
            nEnrDof = nEnrDof * this.getNumberOfDofPerDiscontinuity();
        end

        %------------------------------------------------------------------
        % Obtain the number of discontinuities
        function n = getNumberOfDiscontinuities(this)
            n = size(this.discontinuity,1);
        end

        %------------------------------------------------------------------
        % Calculates the number of degrees of freedom per discontinuity for
        % the enriched element
        function n = getNumberOfDisplacementDofPerDiscontinuity(this)
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
        % Calculates the number of degrees of freedom per discontinuity for
        % the enriched element
        function n = getNumberOfDofPerDiscontinuity(this)
            n = 2;  
        end

        %------------------------------------------------------------------
        % Adds a discontinuity segment to the element
        function addDiscontinuitySegment(this,dseg)
            this.discontinuity = [this.discontinuity; dseg];
        end

        %------------------------------------------------------------------
        % Function to get the displacement jump dofs of this element
        function a = getDisplacementJumpDofs(this,ndofa)
            i = this.nglu + this.nglp + 1;
            j = this.nglu + this.nglp + ndofa;
            a = this.ue(i:j);
        end

        %------------------------------------------------------------------
        % Function to get the displacement jump dofs of this element
        function pj = getPressureJumpDofs(this,ndofa)
            i = this.nglu + this.nglp + ndofa + 1;
            pj = this.ue(i:end);
        end

        %------------------------------------------------------------------
        % Adds the discontinuities dofs to the element dof vector
        function addEnrichmentToDofVector(this)
            nDiscontinuities = this.getNumberOfDiscontinuities();
            dof_a = [];
            for i = 1:nDiscontinuities
                dof_a = [dof_a, this.discontinuity(i).dof];
            end
            this.gle = [this.gle, dof_a];
            this.ngle = length(this.gle);
        end

        %------------------------------------------------------------------
        % Compute the enrichment shape function matrix
        function Nenr = enrichedShapeFncMtrx(this, N, Xcar) 
            nDiscontinuities  = this.getNumberOfDiscontinuities();
            Nenr = zeros(1,nDiscontinuities);  % Each discontinuity has a pressure jump dof
            for i = 1:nDiscontinuities  
                Nenr_i = this.enrichedShapeFncValues(i, N, Xcar);
                Nenr(i) = Nenr_i(1);
            end
        end

        %------------------------------------------------------------------
        % Compute the enrichment shape function values of a discontinuity
        function Nenr_i = enrichedShapeFncValues(this, id, N, Xcar) 
            phi = 0.0;
            h = this.discontinuity(id).heaviside(Xcar);
            for j = 1:this.nnd_el
                Xj = this.node(j,:);
                hj = this.discontinuity(id).heaviside(Xj);
                if (hj > 0.0)
                    phi = phi + N(:, j);
                end
            end
            Nenr_i = [h-phi;
                       -phi;        % Nbot
                      1.0-phi];     % Ntop
        end

        %------------------------------------------------------------------
        % Compute the gradient enrichment matrix
        function G = Gmatrix(this, Bu) 
            nDiscontinuities  = this.getNumberOfDiscontinuities();
            G = zeros(2,nDiscontinuities);  % Each discontinuity has a pressure jump dof
            for i = 1:nDiscontinuities    
                for j = 1:this.nnd_el
                    Xj = this.node(j,:);
                    h = this.discontinuity(i).heaviside(Xj);
                    if (h > 0.0)
                        G(:,i) = G(:,i) - Bu(:, j);
                    end
                end
            end
        end

        %------------------------------------------------------------------
        % Computed the kinematic enrichment matrix for an enriched finite
        % element.
        % It calculated the enrichment matrix by considering the
        % contributions of discontinuities in the element. The enrichment
        % includes translation, stretching and relative rotation modes
        % depending on the configuration
        function Gc = kinematicEnrichment(this, Bu) 
            nDiscontinuities  = this.getNumberOfDiscontinuities();
            nDofDiscontinuity = this.getNumberOfDisplacementDofPerDiscontinuity();
            Gc = zeros(4,nDofDiscontinuity * nDiscontinuities);
            for i = 1:nDiscontinuities    
                Gci = zeros(4,nDofDiscontinuity);
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
                Gc(:,cols) = Gci;
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
        % Compute the enrichment auxiliary function
        function Gv = equilibriumEnrichment(this, Xn)
            
            % Initialize variables
            nDiscontinuities  = this.getNumberOfDiscontinuities();
            nDofDiscontinuity = this.getNumberOfDisplacementDofPerDiscontinuity();
            jumpOrder         = this.displacementJumpOrder();
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
                Gv(:,cols) = Gi;

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
