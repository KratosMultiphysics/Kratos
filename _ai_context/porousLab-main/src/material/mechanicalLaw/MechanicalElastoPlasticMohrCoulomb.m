%% MechanicalElastoPlasticMohrCoulomb Class
% Mohr–Coulomb elastoplastic law (smooth p–q form) with nonassociated flow.
% Uses your base MechanicalElastoPlastic return mapping.
%
% Yield: F = q + (sin(phi)/3)*I1 - c*cos(phi)
% Flow : g = q + (sin(psi)/3)*I1
%
% Depends on PorousMedia fields:
%   material.cohesion        [Pa]
%   material.frictionAngle   [rad]  (phi)
%   material.dilationAngle   [rad]  (psi)
%   material.Kp              [Pa]   (isotropic plastic modulus; 0 => perfect plasticity)
%
% Voigt order: [xx yy zz xy] in plane problems (your IntPoint sets nVar=4).
% The base class handles plane-stress/plane-strain via ip.anm and Ce/De.
%
% Reference:
% de Souza Neto, Eduardo A., Djordje Peric, and David RJ Owen. 
% Computational methods for plasticity: theory and applications.
% John Wiley & Sons, 2011.

classdef MechanicalElastoPlasticMohrCoulomb < MechanicalElastoPlastic

    methods
        %------------------------------------------------------------------
        function this = MechanicalElastoPlasticMohrCoulomb()
            this = this@MechanicalElastoPlastic();
            this.nstVar = 0;
        end
    end

    %% Model-specific pieces required by MechanicalElastoPlastic
    methods
        %------------------------------------------------------------------
        % Compute the stress vector and the constitutive matrix
        function [stress,Dt] = eval(this,material,ip)
            if strcmp(material.stressIntAlgorithm,'implicit')
                 [stress,Dt] = eval@MechanicalElastoPlastic(this,material,ip);
            elseif strcmp(material.stressIntAlgorithm,'alternative')
                [stress,Dt] = this.alternativeStressIntegration(material,ip);
            else
                disp('Error: the given stress integration algorithm is not available');
                disp('Tags of the methods available: ''implicit'', ''alternative''');
                error('Error: stressIntAlgorithm is not available');
            end
        end

        %------------------------------------------------------------------
        % Compute the stress vector and the constitutive matrix
        function [stress,Dt] = alternativeStressIntegration(this,material,ip)

            % Auxiliar constant
            SMALL = 1e-8;

            % Elastic constants
            K = this.bulkModulus(material);
            G = this.shearModulus(material);

            % MC parameters
            phi  = material.frictionAngle;
            psi  = material.dilationAngle;
            c    = material.cohesion;
            sphi = sin(phi);
            cphi = cos(phi);
            spsi = sin(psi);

            % Volumetric strain
            ev = this.volumetricStrain(ip.strain);

            % Pressure trial
            pressure_trial = K * ev;

            % Deviatoric strains
            s_xx_trial = 2*G*(ip.strain(1) - ev/3) + pressure_trial;
            s_yy_trial = 2*G*(ip.strain(2) - ev/3) + pressure_trial;
            s_zz_trial = 2*G*(ip.strain(3) - ev/3) + pressure_trial;
            s_xy_trial = G*ip.strain(4);

            % Stresses spectral decomposition
            [eig12,EIGPRJ] = this.spectralDecomposition([s_xx_trial,s_yy_trial,s_xy_trial]);

            % Constitutive matrix
            De = this.elasticConstitutiveMatrix(material,ip);

            % Principal trial stresses
            pstrs = [eig12(1); eig12(2); s_zz_trial];

            % Identify max/min/mid principal indices
            [p1, ii] = max(pstrs);
            [p3, jj] = min(pstrs);
            mm = setdiff(1:3,[ii,jj]);
            p2 = pstrs(mm);

            % Trial yield on main plane: Φ = (σ1−σ3) + (σ1+σ3) sinφ − 2 c cosφ
            smc = (p1 - p3) + (p1 + p3)*sphi;
            f =  smc - 2*c*cphi;

            % Check if it is in the elastic state
            if f <= this.returnYieldConditionTol
                Dt = De;
                stress = [s_xx_trial;s_yy_trial;s_zz_trial;s_xy_trial];
                return;
            end

            % ==== PLASTIC STEP ======
            
            % One-vector (main plane) closed form (Box 8.5 with H=0)
            a = 4.0*G*(1.0 + spsi*sphi/3.0) + 4.0*K*spsi*sphi;
            lambda = f / a;

            % Update principal stresses
            s1 = p1 - (2.0 * G * (1.0 + spsi / 3.0) + 2.0 * K * spsi) * lambda;
            s2 = p2 + (4.0 / 3.0 * G - 2.0 * K) * spsi * lambda;
            s3 = p3 + (2.0 * G * (1.0 - spsi / 3.0) - 2.0 * K * spsi) * lambda;

            % Check validity of main plane return
            DELTA = max([abs(s1),abs(s2),abs(s3)])*SMALL;
            if s1+DELTA >= s2 && s2+DELTA >= s3
                pstrs(ii) = s1;
                pstrs(jj) = s3;
                pstrs(mm) = s2;
                stress = EIGPRJ * pstrs;
                Dt = this.consistentTangentMainPlane(K,G,sphi,spsi,stress,ip.strain);
                return
            end

            % Edge selector (Eq. 8.90): S = (1−sinψ)σ1 − 2σ2 + (1+sinψ)σ3
            S = (1 - spsi)*p1 - 2*p2 + (1 + spsi)*p3;
            right = (S >= 0);

            % Two-vector (edge) closed form (Box 8.6 with H=0)
            sigma_a  = (p1 - p3) + (p1 + p3)*sphi;
            if right
                sigma_b = (p1 - p2) + (p1 + p2)*sphi;
                b = 2.0*G*(1 + sphi + spsi - sphi*spsi/3) + 4*K*sphi*spsi;
            else
                sigma_b = (p2 - p3) + (p2 + p3)*sphi;
                b = 2.0*G*(1 - sphi - spsi - sphi*spsi/3) + 4*K*sphi*spsi;
            end
            RHS = [sigma_a - 2*c*cphi;  sigma_b - 2*c*cphi];
            A   = [a b; b a];
            sol = A \ RHS;
            lambda_a = sol(1); lambda_b = sol(2);

            % Update principals on the chosen edge (Box 8.6)
            AUX1 = 2.0*G*(1 + spsi/3) + 2.0*K*spsi;
            AUX2 = (4*G/3 - 2*K)*spsi;
            AUX3 = 2.0*G*(1 - spsi/3) - 2.0*K*spsi;
            if right
                s1 = p1 - AUX1*(lambda_a + lambda_b);
                s2 = p2 + AUX2*lambda_a + AUX3*lambda_b;
                s3 = p3 + AUX3*lambda_a + AUX2*lambda_b;
            else
                s1 = p1 - AUX1*lambda_a + AUX2*lambda_b;
                s2 = p2 + AUX2*lambda_a - AUX1*lambda_b;
                s3 = p3 + AUX3*(lambda_a + lambda_b);
            end

            % Check validity
            DELTA = max([abs(s1),abs(s2),abs(s3)])*SMALL;
            if s1+DELTA >= s2 && s2+DELTA >= s3
                % Valid 2-vector return
                pstrs(ii) = s1;
                pstrs(jj) = s3;
                pstrs(mm) = s2;
                stress = EIGPRJ * pstrs;
                Dt = this.consistentTangentEdge(K,G,sphi,spsi,stress,ip.strain,right);
                return
            end

            % Apex return (Box 8.7 with H=0) — closed form
            if abs(sphi) < eps || abs(spsi) < eps
                error('SUMC_PP apex: invalid φ or ψ = 0.');
            end
            P    = (cphi/sphi)*c;                 % P = cot(phi) * c
            s1 = P; s2 = P; s3 = P;
            pstrs(ii) = s1;
            pstrs(jj) = s3;
            pstrs(mm) = s2;
            stress = EIGPRJ * pstrs;
            Dt = this.consistentTangentApex(EIGPRJ,0.0e-8);
        end

        %------------------------------------------------------------------
        % Consistent tangent matrix for the return at the main plane
        function Dt = consistentTangentMainPlane(this,K,G,sphi,spsi,stress,strain)

            % Spectral decomposition of the strain tensor
            [eig12,EIGPRJ,repeat] = this.spectralDecomposition([strain(1),strain(2),0.5*strain(4)]);
            pstrain = [eig12; strain(3)];

            % Total stresses aligned with the principal strains
            pstrs(1) = EIGPRJ(1,1) * stress(1) + EIGPRJ(2,1) * stress(2) + 2*EIGPRJ(4,1) * stress(4);
            pstrs(2) = EIGPRJ(1,2) * stress(1) + EIGPRJ(2,2) * stress(2) + 2*EIGPRJ(4,2) * stress(4);
            pstrs(3) = stress(3);

            % Identify the directions of maximum and minimum principal trial stresses
            ii = 1; jj = 1;
            pstmax = pstrain(ii);
            pstmin = pstrain(jj);
            for i = 2:3
                if (pstrain(i) >= pstmax)
                    ii = i;
                    pstmax = pstrain(i);
                end
                if (pstrain(i) < pstmin)
                    jj = i;
                    pstmin = pstrain(i);
                end
            end
            if (ii ~= 1 && jj ~= 1), mm = 1;
            elseif (ii ~= 2 && jj ~= 2), mm = 2;
            else, mm = 3;
            end

            CONSTA = 4*G*(1 + (sphi*spsi)/3) + 4*K*(sphi*spsi);
            DENOM  = -CONSTA;
        
            B1 = (2*G*(1 + spsi/3) + 2*K*spsi)/DENOM;
            B2 = ((4*G/3 - 2*K)*spsi)/DENOM;
            B3 = (2*G*(1 - spsi/3) - 2*K*spsi)/DENOM;
        
            DP = zeros(3);
        
            % --- B1 block (row II) ---
            DP(ii,ii) = 2*G*(2/3 + B1*(1 + sphi/3)) + K*(1 + 2*B1*sphi);
            DP(ii,mm) = (1/3)*(3*K - 2*G) * (1 + 2*B1*sphi);
            DP(ii,jj) = 2*G*(-1/3 - B1*(1 - sphi/3)) + K*(1 + 2*B1*sphi);
            
            % --- B2 block (row MM) ---
            DP(mm,ii) = 2*G*(-1/3 - B2*(1 + sphi/3)) + K*(1 - 2*B2*sphi);
            DP(mm,mm) = (4*G/3)*(1 + B2*sphi) + K*(1 - 2*B2*sphi);
            DP(mm,jj) = 2*G*(-1/3 + B2*(1 - sphi/3)) + K*(1 - 2*B2*sphi);
            
            % --- B3 block (row JJ) ---
            DP(jj,ii) = 2*G*(-1/3 - B3*(1 + sphi/3)) + K*(1 - 2*B3*sphi);
            DP(jj,mm) = (1/3)*(3*K - 2*G) * (1 - 2*B3*sphi);
            DP(jj,jj) = 2*G*(2/3 + B3*(1 - sphi/3)) + K*(1 - 2*B3*sphi);

            Dt = this.assembleDt(DP, EIGPRJ, pstrain(1:2), pstrs(1:2), true, repeat);

        end

        %------------------------------------------------------------------
        % Consistent tangent matrix for the return at the main plane
        function Dt = consistentTangentEdge(this,K,G,sphi,spsi,stress,strain,right)

            % Spectral decomposition of the strain tensor
            [eig12,EIGPRJ,repeat] = this.spectralDecomposition([strain(1),strain(2),0.5*strain(4)]);
            pstrain = [eig12; strain(3)];

            % Total stresses aligned with the principal strains
            pstrs(1) = EIGPRJ(1,1) * stress(1) + EIGPRJ(2,1) * stress(2) + 2*EIGPRJ(4,1) * stress(4);
            pstrs(2) = EIGPRJ(1,2) * stress(1) + EIGPRJ(2,2) * stress(2) + 2*EIGPRJ(4,2) * stress(4);
            pstrs(3) = stress(3);

            % Identify the directions of maximum and minimum principal trial stresses
            ii = 1; jj = 1;
            pstmax = pstrain(ii);
            pstmin = pstrain(jj);
            for i = 2:3
                if (pstrain(i) >= pstmax)
                    ii = i;
                    pstmax = pstrain(i);
                end
                if (pstrain(i) < pstmin)
                    jj = i;
                    pstmin = pstrain(i);
                end
            end
            if (ii ~= 1 && jj ~= 1), mm = 1;
            elseif (ii ~= 2 && jj ~= 2), mm = 2;
            else, mm = 3;
            end

            SP = sphi*spsi;
          CONSTA = 4*G*(1 + SP/3) + 4*K*SP;
          if right
            CONSTB = 2*G*(1 + sphi + spsi - SP/3) + 4*K*SP;
          else
            CONSTB = 2*G*(1 - sphi - spsi - SP/3) + 4*K*SP;
          end
        
          DRVAA=-CONSTA; DRVBB=-CONSTA; DRVAB=-CONSTB; DRVBA=-CONSTB;
          invDet = 1/(DRVAA*DRVBB - DRVAB*DRVBA);
        
          AUX1 = 2*G*(1 + spsi/3) + 2*K*spsi;
          AUX2 = (4*G/3 - 2*K)*spsi;
          AUX3 = 2*G*(1 - spsi/3) - 2*K*spsi;
        
          DP = zeros(3);
        
          if right
            % ---- CTMC: rows II, MM, JJ with columns (II, MM, JJ) ----
            DP(ii,ii) = K + (4*G/3) + AUX1*(-DRVAB+DRVBB+DRVAA-DRVBA) ...
                        *(2*G + (2*K + 2*G/3)*sphi)*invDet;
            DP(ii,mm) = K - (2*G/3) + AUX1*( 2*G*(DRVAB-DRVAA) ...
                      + ((-DRVAB+DRVBB+DRVAA-DRVBA)*(2*K+2*G/3) + (DRVBA-DRVBB)*2*G)*sphi )*invDet;
            DP(ii,jj) = K - (2*G/3) + AUX1*( 2*G*(DRVBA-DRVBB) ...
                      + ((-DRVAB+DRVBB+DRVAA-DRVBA)*(2*K+2*G/3) + (DRVAB-DRVAA)*2*G)*sphi )*invDet;
        
            DP(mm,ii) = K - (2*G/3) + (AUX2*(DRVAB-DRVBB) + AUX3*(DRVBA-DRVAA)) ...
                      *(2*G + (2*K + 2*G/3)*sphi)*invDet;
            DP(mm,mm) = K + (4*G/3) + ( AUX2*( (2*K*(DRVAB-DRVBB) + (DRVAB*(2*G/3)+DRVBB*(4*G/3)))*sphi - DRVAB*2*G ) ...
                      + AUX3*( DRVAA*2*G + (2*K*(DRVBA-DRVAA) - (DRVAA*(2*G/3)+DRVBA*(4*G/3)))*sphi ) )*invDet;
            DP(mm,jj) = K - (2*G/3) + ( AUX2*( (2*K*(DRVAB-DRVBB) - (DRVBB*(2*G/3)+DRVAB*(4*G/3)))*sphi + DRVBB*2*G ) ...
                      + AUX3*( (2*K*(DRVBA-DRVAA) + (DRVAA*(4*G/3)+DRVBA*(2*G/3)))*sphi - DRVBA*2*G ) )*invDet;
        
            DP(jj,ii) = K - (2*G/3) + ((AUX2*(DRVBA-DRVAA)+AUX3*(DRVAB-DRVBB)) * ((2*K+2*G/3)*sphi + 2*G))*invDet;
            DP(jj,mm) = K - (2*G/3) + ( AUX2*( (2*K*(DRVBA-DRVAA) - (DRVBA*(4*G/3)+DRVAA*(2*G/3)))*sphi + DRVAA*2*G ) ...
                      + AUX3*( (2*K*(DRVAB-DRVBB) + (DRVAB*(2*G/3)+DRVBB*(4*G/3)))*sphi - DRVAB*2*G ) )*invDet;
            DP(jj,jj) = K + (4*G/3) + ( AUX2*( (2*K*(DRVBA-DRVAA) + (DRVAA*(4*G/3)+DRVBA*(2*G/3)))*sphi - DRVBA*2*G ) ...
                      + AUX3*( (2*K*(DRVAB-DRVBB) - (DRVAB*(4*G/3)+DRVBB*(2*G/3)))*sphi + DRVBB*2*G ) )*invDet;
        
          else
            % ---- CTMC: rows II, MM, JJ with columns (II, MM, JJ) ----
            DP(ii,ii) = K + (4*G/3) + ( AUX1*( (2*K*(DRVBB-DRVAB) + (DRVAB*(4*G/3)+DRVBB*(2*G/3)) )*sphi + DRVBB*2*G ) ...
                      + AUX2*( (2*K*(DRVBA-DRVAA) + (DRVAA*(4*G/3)+DRVBA*(2*G/3)) )*sphi + DRVBA*2*G ) )*invDet;
            DP(ii,mm) = K - (2*G/3) + ( AUX1*( (2*K*(DRVBB-DRVAB) - (DRVAB*(2*G/3)+DRVBB*(4*G/3)) )*sphi - DRVAB*2*G ) ...
                      + AUX2*( (2*K*(DRVBA-DRVAA) - (DRVAA*(2*G/3)+DRVBA*(4*G/3)) )*sphi - DRVAA*2*G ) )*invDet;
            DP(ii,jj) = K - (2*G/3) + ((AUX1*(DRVBB-DRVAB)+AUX2*(DRVBA-DRVAA)) * ((2*K+2*G/3)*sphi - 2*G))*invDet;
        
            DP(mm,ii) = K - (2*G/3) + ( AUX1*( (2*K*(DRVAA-DRVBA) - (DRVAA*(4*G/3)+DRVBA*(2*G/3)) )*sphi - DRVBA*2*G ) ...
                      + AUX2*( (2*K*(DRVAB-DRVBB) - (DRVAB*(4*G/3)+DRVBB*(2*G/3)) )*sphi - DRVBB*2*G ) )*invDet;
            DP(mm,mm) = K + (4*G/3) + ( AUX1*( (2*K*(DRVAA-DRVBA) + (DRVAA*(2*G/3)+DRVBA*(4*G/3)) )*sphi + DRVAA*2*G ) ...
                      + AUX2*( (2*K*(DRVAB-DRVBB) + (DRVAB*(2*G/3)+DRVBB*(4*G/3)) )*sphi + DRVAB*2*G ) )*invDet;
            DP(mm,jj) = K - (2*G/3) + ((AUX1*(DRVAA-DRVBA)+AUX2*(DRVAB-DRVBB)) * ((2*K+2*G/3)*sphi - 2*G))*invDet;
        
            DP(jj,ii) = K - (2*G/3) + AUX3*( (2*K*(DRVAB-DRVBB-DRVAA+DRVBA) + (DRVAA-DRVAB)*(4*G/3) + (DRVBA-DRVBB)*(2*G/3) )*sphi + (DRVBA-DRVBB)*2*G )*invDet;
            DP(jj,mm) = K - (2*G/3) + AUX3*( (2*K*(DRVAB-DRVBB-DRVAA+DRVBA) + (DRVAB-DRVAA)*(2*G/3) + (DRVBB-DRVBA)*(4*G/3) )*sphi + (DRVAB-DRVAA)*2*G )*invDet;
            DP(jj,jj) = K + (4*G/3) + AUX3*( (DRVAB-DRVBB-DRVAA+DRVBA) * ((2*K+2*G/3)*sphi - 2*G) )*invDet;
          end
        
          Dt = this.assembleDt(DP, EIGPRJ, pstrain(1:2), pstrs(1:2), true, repeat);
        
        end

         %------------------------------------------------------------------
        function DYDX = assembleDt(~,DEIGY, EIGPRJ, eigX, eigY, OUTOFP, REPEAT)
        % DGISO2 (Souza Neto) adapted to Voigt order [xx yy zz xy]
        % In-plane indices = [1 2 4] (xx, yy, xy); out-of-plane = 3 (zz).
        % Inputs:
        %   DEIGY  [3x3]  dY_principal / dX_principal in principal basis
        %   EIGPRJ [4x2]  projector columns for the two in-plane principal dirs
        %                 (rows in [xx yy zz xy] order; row 3 (zz) must be 0)
        %   eigX   [2x1]  in-plane principal values of X
        %   eigY   [2x1]  in-plane principal values of Y
        %   OUTOFP logical  include out-of-plane (zz) coupling
        %   REPEAT logical  repeated (or nearly) in-plane eigenvalues
        %
        % Output:
        %   DYDX   [4x4]  derivative dY/dX in Voigt [xx yy zz xy]
        
            % index maps (Voigt: [xx yy zz xy])
            ip  = [1 2 4];   % in-plane (xx, yy, xy)
            oop = 3;         % out-of-plane (zz)
        
            % in-plane identity (engineering shear => 0.5 on (xy,xy))
            FOID3 = diag([1, 1, 0.5]);   % 3x3 acting on [xx yy xy]
        
            % "sum" projector used by formulas (1 on xx,yy; 0 on zz,xy)
            SOPID = [1; 1; 0; 0];        % 4x1 in [xx yy zz xy]
        
            % selector for out-of-plane component (zz)
            EIGPR3 = zeros(4,1); EIGPR3(oop) = 1;
        
            DYDX = zeros(4,4);
        
            if REPEAT
                % -------- repeated in-plane eigenvalues ----------
                % in-plane 3x3 block (map into {1,2,4})
                for I = 1:3
                    for J = 1:3
                        DYDX(ip(I), ip(J)) = ...
                            (DEIGY(1,1) - DEIGY(1,2))*FOID3(I,J) + ...
                             DEIGY(1,2) * (SOPID(ip(I)) * SOPID(ip(J)));
                    end
                end
        
                if OUTOFP
                    % terms touching zz (oop)
                    for i = 1:4
                        for j = 1:4
                            if (i == oop) || (j == oop)
                                DYDX(i,j) = ...
                                      DEIGY(1,3) * SOPID(i)  * EIGPR3(j) ...
                                    + DEIGY(3,1) * EIGPR3(i) * SOPID(j)  ...
                                    + DEIGY(3,3) * EIGPR3(i) * EIGPR3(j);
                            end
                        end
                    end
                end
        
            else
                % -------- distinct in-plane eigenvalues ----------
                A1 = (eigY(1) - eigY(2)) / (eigX(1) - eigX(2));
        
                % in-plane 3x3 block (map into {1,2,4})
                for I = 1:3
                    for J = 1:3
                        e11 = EIGPRJ(ip(I),1); e21 = EIGPRJ(ip(J),1);
                        e12 = EIGPRJ(ip(I),2); e22 = EIGPRJ(ip(J),2);
                        DYDX(ip(I), ip(J)) = ...
                            A1 * ( FOID3(I,J) - e11*e21 - e12*e22 ) ...
                          + DEIGY(1,1)*e11*e21 + DEIGY(1,2)*e11*e22 ...
                          + DEIGY(2,1)*e12*e21 + DEIGY(2,2)*e12*e22;
                    end
                end
        
                if OUTOFP
                    % coupling with zz (oop)
                    for i = 1:4
                        for j = 1:4
                            if (i == oop) || (j == oop)
                                DYDX(i,j) = ...
                                      DEIGY(1,3) * EIGPRJ(i,1) * EIGPR3(j) ...
                                    + DEIGY(2,3) * EIGPRJ(i,2) * EIGPR3(j) ...
                                    + DEIGY(3,1) * EIGPR3(i)   * EIGPRJ(j,1) ...
                                    + DEIGY(3,2) * EIGPR3(i)   * EIGPRJ(j,2) ...
                                    + DEIGY(3,3) * EIGPR3(i)   * EIGPR3(j);
                            end
                        end
                    end
                end
            end
        end

        function [eig12,EIGPRJ,repeatFlag] = spectralDecomposition(~,X)
            % eigenvalues (closed form)
            trx = X(1) + X(2);
            B   = sqrt((X(1) - X(2))^2 + 4*X(3)*X(3));
            eig1 = 0.5*(trx + B);
            eig2 = 0.5*(trx - B);
            eig12 = [eig1; eig2];
            
            % degeneracy check (same normalization as SPDEC2)
            SMALL = 1e-5;
            diffn = abs(eig1 - eig2);
            amax  = max(abs(eig1), abs(eig2));
            if amax ~= 0, diffn = diffn/amax; end
            repeatFlag = (diffn < SMALL);
            
            EIGPRJ = zeros(4,3);
            if repeatFlag
                % robust fallback: compute eigenvectors numerically (like JACOB)
                A = [X(1), X(3); X(3), X(2)];
                [V, D] = eig((A + A.')/2);
                % sort by descending eigenvalue
                [dvals, idx] = sort(diag(D), 'descend'); %#ok<ASGLU>
                V = V(:, idx);
                v1 = V(:,1); v2 = V(:,2);
            
                % projectors v v^T in Voigt [xx; yy; zz; xy], with zz row=0
                EIGPRJ(:,1) = [v1(1)^2; v1(2)^2; 0; v1(1)*v1(2)];
                EIGPRJ(:,2) = [v2(1)^2; v2(2)^2; 0; v2(1)*v2(2)];
            else
                % closed-form projectors (exactly SPDEC2’s formula)
                for k = 1:2
                    lam = eig12(k);
                    b   = lam - trx;                 % = - other eigenvalue
                    c   = 1.0 / (lam + b);           % = 1/(lam - other)
                    EIGPRJ(1,k) = c*(X(1) + b);
                    EIGPRJ(2,k) = c*(X(2) + b);
                    EIGPRJ(4,k) = c*(X(3));
                    EIGPRJ(3,k) = 0.0;               % zz component not in 2D projector
                end
            end
            EIGPRJ(3,3) = 1;
        end


        %------------------------------------------------------------------
        % Consistent tangent matrix for the return at the main plane
        function Dt = consistentTangentApex(~,EIGPRJ,reg)
            if nargin < 3, reg = 0; end     % optional tiny regularization
            DP = reg*eye(3);                 % H=0 ⇒ true DP would be zeros(3)
            Dt = EIGPRJ * DP * EIGPRJ.';
        end

        %------------------------------------------------------------------
        function Dt = rotateDPtoCartesian(~, DP, strain, sprinc)
            % --- projectors from TRIAL STRAIN (Souza-Neto CTMC) ---
            exx = strain(1); eyy = strain(2); exy = strain(4)/2;   % tensor shear
            A = [exx, exy; exy, eyy];
            [V,D] = eig((A+A.')/2);
            [~,idx] = sort(diag(D),'descend');
            V = V(:,idx);
            v1 = V(:,1); v2 = V(:,2);
            e1 = D(idx(1),idx(1)); e2 = D(idx(2),idx(2));
        
            % --- 4x3 projectors: stress-side (sP) and strain-side (eP) ---
            % Voigt order [xx yy zz xy]; stress uses τ_xy; strain uses γ_xy = 2ε_xy
            sP = zeros(4,3);   eP = zeros(4,3);
        
            % σ1/ε1 direction
            sP(:,1) = [v1(1)^2; v1(2)^2; 0; v1(1)*v1(2)];
            eP(:,1) = [v1(1)^2; v1(2)^2; 0; 2*v1(1)*v1(2)];
        
            % σ2/ε2 direction
            sP(:,2) = [v2(1)^2; v2(2)^2; 0; v2(1)*v2(2)];
            eP(:,2) = [v2(1)^2; v2(2)^2; 0; 2*v2(1)*v2(2)];
        
            % out-of-plane (zz)
            sP(:,3) = [0;0;1;0];
            eP(:,3) = [0;0;1;0];
        
            % --- base push-forward ---
            % δσ = sP * DP * (eP^T δε)
            Dt = sP * DP * eP.';
        
            % --- projector-derivative (rotation) term (DGISO2) ---
            denom = max(abs(e1 - e2), 1e-12);
            alpha = (sprinc(1) - sprinc(2)) / denom;
        
            % note: stress-side projector on left, strain-side projector on right
            Crot = alpha * ( sP(:,1)*eP(:,2).' + sP(:,2)*eP(:,1).' );
            Dt = Dt + Crot;
        end

        %------------------------------------------------------------------
        % Yield function F(sigma)
        function f = yieldCondition(this, material, ~, stress, ~)

            J2 = max(this.stressInvariantJ2(stress), 1e-10);
            J3 = this.stressInvariantJ3(stress);
            I1 = this.stressInvariantI1(stress);
            J  = sqrt(J2);
        
            x     = -(3.0*sqrt(3.0)/2.0)*(J3/(J^3));             % sin(3θ)
            x     = max(min(x, 1-1e-10), -1+1e-10);
            theta = asin(x)/3.0;
        
            phi = material.frictionAngle;  c = material.cohesion;
            A   = cos(theta) - sin(theta)*sin(phi)/sqrt(3.0);
        
            f = (sin(phi)/3.0)*I1 + J*A - c*cos(phi);
        end

        %------------------------------------------------------------------
        % dF/dsigma
        function df = yieldStressGradient(this, material, ~, stress, ~)
            
            % Material parameters
            phi = material.frictionAngle;

            % invariants and gradients
            J2 = this.stressInvariantJ2(stress);
            J3 = this.stressInvariantJ3(stress);
            dI1 = this.gradientI1();
            dJ2 = this.gradientJ2(stress);
            dJ3 = this.gradientJ3(stress);
        
            % robust scalars
            epsJ2 = 1.0e-10;
            J2v = max(J2, epsJ2);
            J    = sqrt(J2v);
        
            % angle machinery reusing x and cos(3θ)
            x = -(3.0*sqrt(3.0)/2.0) * (J3/(J^3));       % sin(3θ)
            x = max(min(x, 1 - 1.0e-10), -1 + 1.0e-10);
            theta = asin(x)/3.0;
            cos3  = sqrt(max(0.0, 1.0 - x*x));          % = cos(3θ) on principal branch

            % Auxiliary values
            A = cos(theta) - sin(theta) * sin(phi) / sqrt(3.0);
            dA = -sin(theta) - cos(theta) * sin(phi) / sqrt(3.0);

            % Derivatives of the yield function wrt the invariants
            dfdI1 = sin(phi)/3.0;
            % tan(3θ) = x / cos(3θ)
            if cos3 > 1.0e-12
                dfdJ2 = (A - dA*(x/cos3)) / (2.0*J);
                dfdJ3 = -sqrt(3.0)*dA / (2.0*J2v*cos3);
            else
                % freeze θ-dependence at the cap
                dfdJ2 = A / (2.0*J);
                dfdJ3 = 0.0;
            end

            % Gradient
            df = dfdI1 * dI1 + dfdJ2 * dJ2 + dfdJ3 * dJ3;
        end

        %------------------------------------------------------------------
        % Flow vector n = ∂g/∂σ  (nonassociated if psi ≠ phi)
        function n = flowVector(this, material, ~, stress, ~)
            
            % Material parameters
            psi = material.dilationAngle;

            % invariants and gradients
            J2 = this.stressInvariantJ2(stress);
            J3 = this.stressInvariantJ3(stress);
            dI1 = this.gradientI1();
            dJ2 = this.gradientJ2(stress);
            dJ3 = this.gradientJ3(stress);
        
            % robust scalars
            epsJ2 = 1.0e-10;
            J2v = max(J2, epsJ2);
            J    = sqrt(J2v);
        
            % angle machinery reusing x and cos(3θ)
            x = -(3.0*sqrt(3.0)/2.0) * (J3/(J^3));       % sin(3θ)
            x = max(min(x, 1 - 1.0e-10), -1 + 1.0e-10);
            theta = asin(x)/3.0;
            cos3  = sqrt(max(0.0, 1.0 - x*x));          % = cos(3θ) on principal branch

            % Auxiliary values
            A = cos(theta) - sin(theta) * sin(psi) / sqrt(3.0);
            dA = -sin(theta) - cos(theta) * sin(psi) / sqrt(3.0);

            % Derivatives of the yield function wrt the invariants
            dfdI1 = sin(psi)/3.0;
            % tan(3θ) = x / cos(3θ)
            if cos3 > 1.0e-12
                dfdJ2 = (A - dA*(x/cos3)) / (2.0*J);
                dfdJ3 = -sqrt(3.0)*dA / (2.0*J2v*cos3);
            else
                % freeze θ-dependence at the cap
                dfdJ2 = A / (2.0*J);
                dfdJ3 = 0.0;
            end

            % Gradient
            n = dfdI1 * dI1 + dfdJ2 * dJ2 + dfdJ3 * dJ3;
        end

        %------------------------------------------------------------------
        % dn/dsigma = ∂^2 g / ∂σ∂σ (Hessian of potential)
        function dn = flowStressGradient(this, material, ~, stress, ~)
            % Material parameters
            psi = material.dilationAngle;

            % Invariants and its gradients
            J2 = this.stressInvariantJ2(stress);
            J3 = this.stressInvariantJ3(stress);
            dJ2 = this.gradientJ2(stress);
            dJ3 = this.gradientJ3(stress);
            ddJ2 = this.hessianJ2();
            ddJ3 = this.hessianJ3(stress);

            % Robust scalars
            epsJ2 = 1.0e-10;  tol = 1.0e-10;  sqrt3 = sqrt(3.0);
            J2v = max(J2, epsJ2);
            J   = sqrt(J2v);
        
            % Angle machinery: x = sin(3θ), clamp once
            x    = -(3.0*sqrt3/2.0) * (J3 / (J^3));         % sin(3θ)
            x    = max(min(x, 1 - tol), -1 + tol);
            theta = asin(x)/3.0;
            cos3  = sqrt(max(0.0, 1.0 - x*x));              % = cos(3θ)
            tan3  = x / max(cos3, 1.0e-16);
        
            % A(θ), A'(θ), A''(θ) with dilation angle ψ
            A   = cos(theta) - sin(theta)*sin(psi)/sqrt3;
            dA  = -sin(theta) - cos(theta)*sin(psi)/sqrt3;
            d2A = -cos(theta) + sin(theta)*sin(psi)/sqrt3;
        
            % Coefficients (∂g/∂J2, ∂g/∂J3, and second partials)
            if cos3 > 1.0e-12
                C2  = (A - dA*tan3) / (2.0*J);                  % ∂g/∂J2
                C3  = -sqrt3 * dA / (2.0 * J2v * cos3);         % ∂g/∂J3
                C4  = d2A + 3.0 * tan3 * dA;                    % helper
                C22 = -(A - tan3*tan3*C4 - 3.0*tan3*dA) / (4.0 * J^3);
                C23 =  (0.5*tan3*C4 + dA) * sqrt3 / (2.0 * J2v^2 * cos3);
                C33 =  3.0 * C4 / (4.0 * J2v^(2.5) * cos3 * cos3);
            else
                % Freeze θ-dependence at clamp for stability
                C2  = A / (2.0*J);
                C3  = 0.0;
                C22 = -A / (4.0*J^3);
                C23 = 0.0;
                C33 = 0.0;
            end
            
            % Flow vector gradient
            dn = C2 * ddJ2 + C3 * ddJ3 ...
               + C22 * (dJ2 * dJ2.') ...
               + C23 * (dJ2 * dJ3.' + dJ3 * dJ2.') ...
               + C33 * (dJ3 * dJ3.');
        end

        %------------------------------------------------------------------
        % Gradient of the yield function wrt to the state variables vector
        function dfda = yieldStateGradient(~,~,~,~,~)
            dfda = zeros(0,1);
        end

        %------------------------------------------------------------------
        % Flow vector gradient wrt to the state variables vector
        function dnda = flowStateGradient(~,~,ip,~,~)
            dnda = zeros(ip.nVar,0);
        end

        %------------------------------------------------------------------
        % Internal/state variables evolution
        function h = stateEvolution(~, ~, ~, ~, ~)
            h = zeros(0,1);
        end

        %------------------------------------------------------------------
        % Gradient of the internal/state variables law wrt to the stress vector
        function dhds = stateStressGradient(~, ~, ip, ~, ~)
            dhds = zeros(0,ip.nVar);
        end

        %------------------------------------------------------------------
        % Gradient of the internal/state variables law wrt to the state variables
        function dhda = stateStateGradient(~, ~, ~, ~, ~)
            dhda = zeros(0,0);
        end
    end

    methods(Static)
        function arr_new = place_back(arr, idx)
            arr_new = zeros(3,1);
            arr_new(idx(1)) = arr(1);   % S1 -> index of max trial principal
            arr_new(idx(2)) = arr(2);   % S2 -> middle
            arr_new(idx(3)) = arr(3);   % S3 -> min
        end
    end
end
