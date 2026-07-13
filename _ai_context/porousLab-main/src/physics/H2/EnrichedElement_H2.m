%% EnrichedElement_H2 class
% This class extends the _RegularElement_H2_ class to define a finite 
% element for single-phase fluid flow that incorporates enriched elements 
% to handle discontinuities. It provides methods to compute element data, 
% manage discontinuities, and calculate enriched degrees of freedom.
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
classdef EnrichedElement_H2 < RegularElement_H2   
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        discontinuity = [];
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = EnrichedElement_H2(node, elem, t, ...
                mat, intOrder, glp, glpg, massLumping, lumpStrategy, ...
                isAxisSymmetric)
            this = this@RegularElement_H2(node, elem, t, ...
                mat, intOrder, glp, glpg, massLumping, lumpStrategy, ...
                isAxisSymmetric);
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
           
           [Ke, Ce, fi, fe, dfidu] = elementData@RegularElement_H2(this);

           if isempty(this.discontinuity)
               return
           else
               [Ke, Ce, fi, fe, dfidu] = this.addDiscontinuityData(Ke, Ce, fi, fe, dfidu);
           end
        end

        %------------------------------------------------------------------
        % Computes and assembles the sub-matrices and sub-vectors for an
        % enriched finite element
        function [Ke, Ce, fi, fe, dfidu] = addDiscontinuityData(this, Ke, Ce, fi, fe, dfidu)
             % Compute contribution of the discontinuity
              [Kd, Cd, fid, fed, dfiddu, T] = getDiscontinuitiesData(this);

              % Add contribution of the discontinuities
              Ke    = Ke    + T' * Kd * T;
              Ce    = Ce    + T' * Cd * T;
              fi    = fi    + T' * fid;
              fe    = fe    + T' * fed;
              dfidu = dfidu + T' * dfiddu * T;
        end

        %------------------------------------------------------------------
        % Computes the stiffness matrix and force vector contributions
        % from discontinuities in the enriched element
        function [Kd, Cd, fid, fed, dfiddu, T] = getDiscontinuitiesData(this)

            nDiscontinuities  = this.getNumberOfDiscontinuities();
            nPIntDofs         = nDiscontinuities * 2;

            % Initialize the output data 
            Kd = zeros(2*nPIntDofs,2*nPIntDofs);
            Cd = zeros(2*nPIntDofs,2*nPIntDofs);
            fid = zeros(2*nPIntDofs,1);
            fed = zeros(2*nPIntDofs,1);
            dfiddu = zeros(2*nPIntDofs,2*nPIntDofs);

            % Condensation matrix
            Td = zeros(nPIntDofs,this.nglp);

            % Vector of the nodal pore-pressure dofs
            pl = this.getNodalLiquidPressure();
            pg = this.getNodalGasPressure();

            % Vector with the old nodal dofs
            plOld = this.getOldNodalLiquidPressure(); 
            pgOld = this.getOldNodalGasPressure();

            % Initialize the discontinuity nodal dofs
            pld    = zeros(2,1);
            pgd    = zeros(2,1);
            pldOld = zeros(2,1);
            pgdOld = zeros(2,1);

            % Loop through the discontinuities
            for i = 1:nDiscontinuities

                % Internal pressure dofs associated with this discontinuity segment
                pint_dofs = 4*(i-1)+1 : 4*i;
                d_dofs = 2*(i-1)+1 : 2*i;

                % Loop through the nodes of the discontinuity to fill the
                % condensation matrix
                for j = 1:2
                    X = this.discontinuity(i).node(j,:);
                    Xn = this.shape.coordCartesianToNatural(this.node,X);
                    Np = this.shape.shapeFncMtrx(Xn);
                    pld(j) = Np * pl;
                    pgd(j) = Np * pg;
                    pldOld(j) = Np * plOld;
                    pgdOld(j) = Np * pgOld;
                    Td(d_dofs(j),:) = Np;
                end

                % Get the discontinuity data
                [Kdi, Cdi, fidi, fedi, dfiddui] = this.discontinuity(i).elementData(pld, pgd, pldOld, pgdOld, this.g, this.gravityOn, this.DTime);

                % Assemble the contribution of this discontinuity
                Kd(pint_dofs,pint_dofs) = Kd(pint_dofs,pint_dofs) + Kdi;
                Cd(pint_dofs,pint_dofs) = Cd(pint_dofs,pint_dofs) + Cdi;
                fid(pint_dofs) = fid(pint_dofs) + fidi;
                fed(pint_dofs) = fed(pint_dofs) + fedi;
                dfiddu(pint_dofs,pint_dofs) = dfiddu(pint_dofs,pint_dofs) + dfiddui;
    
            end

            % Initialize the condensation matrix
            T = [ Td                         , zeros(nPIntDofs, this.nglp);
                  zeros(nPIntDofs, this.nglp), Td];
            
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
        function n = getNumberOfDofPerDiscontinuity(~)
            n = 2;  
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
        
    end
end
