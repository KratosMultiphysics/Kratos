%% MechanicalLaw Class
% This is an abstract class used to implement mechanical constitutive
% laws. It provides methods for evaluating stress, strain, and their
% invariants, as well as utilities for handling plane stress conditions
% and tensor operations.
% 
%% Methods
% * *hydrostaticStress*: Computes the hydrostatic stress.
% * *deviatoricStress*: Computes the deviatoric stress.
% * *stressInvariantI1*: Computes the first stress invariant (I1).
% * *stressInvariantI2*: Computes the second stress invariant (I2).
% * *stressInvariantJ2*: Computes the second deviatoric stress invariant 
%                       (J2).
% * *gradientI1*: Computes the gradient of the first stress invariant (I1).
% * *gradientJ2*: Computes the gradient of the second deviatoric stress 
%                invariant (J2).
% * *hessianJ2*: Computes the Hessian matrix of the second deviatoric 
%               stress invariant (J2).
% * *vonMisesStress*: Computes the von Mises stress.
% * *vonMisesStressGradient*: Computes the gradient of the von Mises stress.
% * *vonMisesStressHessian*: Computes the Hessian matrix of the von Mises 
%                           stress.
% * *planeStressConstitutiveMatrix*: Computes the plane stress constitutive 
%                                   matrix.
% * *elasticOutOfPlaneStrain*: Computes the elastic out-of-plane strain for 
%                             plane stress conditions.
% * *normTensor*: Computes the norm of a strain tensor.
% * *strainInvariantI1*: Computes the first strain invariant (I1).
% * *volumetricStrain*: Computes the volumetric strain.
% * *deviatoricStrain*: Computes the deviatoric strain.
% * *normDeviatoricStrain*: Computes the norm of the deviatoric strain.
% * *strainInvariantI2*: Computes the second strain invariant (I2).
% * *strainInvariantJ2*: Computes the second deviatoric strain invariant 
%                       (J2).
% * *gradientStrainInvariantI1*: Computes the gradient of the first strain 
%                               invariant (I1).
% * *gradientStrainInvariantJ2*: Computes the gradient of the second 
%                               deviatoric strain invariant (J2).
% * *fourthOrderSymTensor*: Returns the fourth-order symmetric tensor.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class Definition
classdef MechanicalLaw < handle  
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        nstVar = 0;   % Number of state variables
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MechanicalLaw()
        end
    end

    %% Abstract methods
    methods(Abstract)

        [stress,De] = eval(this,material,ip)

    end

    %% Public methods
    methods
        %% Stress utilities
        % Considering a 2D stress tensor: Stress = [sxx, syy, szz, sxy]

        %------------------------------------------------------------------
        % Principal stresses
        function [eigx,EIGPRJ] = principalStresses(~,stress)

            sxx = stress(1);
            syy = stress(2);
            sxy = stress(4);

            % eigenvalues (closed form)
            trx = sxx + syy;
            B   = sqrt((sxx - syy)^2 + 4*sxy*sxy);
            eig1 = 0.5*(trx + B);
            eig2 = 0.5*(trx - B);
            eigx = [eig1; eig2];
            
            % degeneracy check (same normalization as SPDEC2)
            SMALL = 1e-5;
            diffn = abs(eig1 - eig2);
            amax  = max(abs(eig1), abs(eig2));
            if amax ~= 0, diffn = diffn/amax; end
            repeatFlag = (diffn < SMALL);
            
            EIGPRJ = zeros(4,3);
            if repeatFlag
                % robust fallback: compute eigenvectors numerically (like JACOB)
                A = [sxx, sxy; sxy, syy];
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
                    lam = eigx(k);
                    b   = lam - trx;                 % = - other eigenvalue
                    c   = 1.0 / (lam + b);           % = 1/(lam - other)
                    EIGPRJ(1,k) = c*(sxx + b);
                    EIGPRJ(2,k) = c*(syy + b);
                    EIGPRJ(4,k) = c*(sxy);
                    EIGPRJ(3,k) = 0.0;               % zz component not in 2D projector
                end
            end
            EIGPRJ(3,3) = 1;
        end

        %------------------------------------------------------------------
        function stress = cartesianStressesFromPrincipal(~, p, Q, szz)
            n1 = Q(:,1);
            n2 = Q(:,2);
            STRESS = p(1)*(n1*n1') + p(2)*(n2*n2');
            stress = [STRESS(1,1);
                      STRESS(2,2);
                      szz;
                      STRESS(1,2)];
        end
        
        %------------------------------------------------------------------
        % Hydrostatic stress
        function sh = hydrostaticStress(this,stress)
            sh = this.stressInvariantI1(stress) / 3.0;
        end

        %------------------------------------------------------------------
        % Deviatoric stress
        function sd = deviatoricStress(this,stress)
            sh = this.hydrostaticStress(stress);
            Im = [1.0;1.0;1.0;0.0];
            sd = stress - sh * Im;
        end

        %------------------------------------------------------------------
        % I1 stress invariant
        function I1 = stressInvariantI1(~,stress)
            I1 = stress(1) + stress(2) + stress(3);
        end

        %------------------------------------------------------------------
        % I2 stress invariant
        function I2 = stressInvariantI2(~,stress)
            sxx = stress(1);
            syy = stress(2);
            szz = stress(3);
            sxy = stress(4);
            I2 =  sxx*syy + syy*szz + sxx*szz - sxy*sxy;
        end

        %------------------------------------------------------------------
        % I3 stress invariant
        function I3 = stressInvariantI3(~,stress)
            sxx = stress(1);
            syy = stress(2);
            szz = stress(3);
            sxy = stress(4);
            I3 =  sxx * syy * szz - szz * sxy * sxy;
        end

        %------------------------------------------------------------------
        % J2 stress invariant
        function J2 = stressInvariantJ2(this,stress)
            I1 = this.stressInvariantI1(stress);
            I2 = this.stressInvariantI2(stress);
            J2  = I1*I1/3.0 - I2;
            J2  = max(J2,0.0);
        end

        %------------------------------------------------------------------
        % J3 stress invariant
        function J3 = stressInvariantJ3(this,stress)
            I1 = this.stressInvariantI1(stress);
            I2 = this.stressInvariantI2(stress);
            I3 = this.stressInvariantI3(stress);
            J3  = 2.0 * I1 * I1 * I1 / 27.0 - I1 * I2 / 3.0 + I3;
        end

        %------------------------------------------------------------------
        % Lode angle
        function theta = lodeAngle(this, stress)
            theta = 0.0;
            J2 = this.stressInvariantJ2(stress);
            if J2 < 1.0E-8
                return
            end
            J3 = this.stressInvariantJ3(stress);
            x = -3.0 * sqrt(3.0) * J3 / (2.0 * J2 * sqrt(J2));
            x = max(min(x, 1 - 1.0E-10), -1 + 1.0E-10);
            theta = asin(x) / 3.0;
        end

        %------------------------------------------------------------------
        % Gradient of the I1 stress invariant
        function dI1 = gradientI1(~,~)
            dI1 = [1.0;1.0;1.0;0.0];
        end

        %------------------------------------------------------------------
        % Gradient of the J2 stress invariant
        function dJ2 = gradientJ2(~,stress)
            sxx = stress(1);
            syy = stress(2);
            szz = stress(3);
            sxy = stress(4);
            dJ2 = zeros(4,1);
            dJ2(1) = (2.0 * sxx - syy - szz)/3.0;
            dJ2(2) = (2.0 * syy - sxx - szz)/3.0;
            dJ2(3) = (2.0 * szz - syy - sxx)/3.0;
            dJ2(4) = 2.0 * sxy;
        end

        %------------------------------------------------------------------
        % Gradient of the J3 stress invariant
        function dJ3 = gradientJ3(~,stress)
            sxx = stress(1);
            syy = stress(2);
            szz = stress(3);
            sxy = stress(4);
            dJ3 = zeros(4,1);
            dJ3(1) = (2*sxx^2)/9 - (2*sxx*syy)/9 - (2*sxx*szz)/9 + sxy^2/3 - syy^2/9 + (4*syy*szz)/9 - szz^2/9;
            dJ3(2) =  - sxx^2/9 - (2*sxx*syy)/9 + (4*sxx*szz)/9 + sxy^2/3 + (2*syy^2)/9 - (2*syy*szz)/9 - szz^2/9;
            dJ3(3) =  - sxx^2/9 + (4*sxx*syy)/9 - (2*sxx*szz)/9 - (2*sxy^2)/3 - syy^2/9 - (2*syy*szz)/9 + (2*szz^2)/9;
            dJ3(4) = (2*sxy*(sxx + syy - 2*szz))/3;
        end

        %------------------------------------------------------------------
        % Hessian of the J2 stress invariant
        function d2J2 = hessianJ2(~)
            d2J2 = [  2.0 , -1.0 , -1.0 , 0.0;
                     -1.0 ,  2.0 , -1.0 , 0.0;
                     -1.0 , -1.0 ,  2.0 , 0.0;
                      0.0 ,  0.0  , 0.0 , 6.0 ]/3.0;
        end

        %------------------------------------------------------------------
        % Hessian of the J3 stress invariant
        function d2J3 = hessianJ3(~,stress)
            sxx = stress(1);
            syy = stress(2);
            szz = stress(3);
            sxy = stress(4);
            d2J3 = [(4*sxx)/9 - (2*syy)/9 - (2*szz)/9, (4*szz)/9 - (2*syy)/9 - (2*sxx)/9, (4*syy)/9 - (2*sxx)/9 - (2*szz)/9, (2*sxy)/3;
                    (4*szz)/9 - (2*syy)/9 - (2*sxx)/9, (4*syy)/9 - (2*sxx)/9 - (2*szz)/9, (4*sxx)/9 - (2*syy)/9 - (2*szz)/9, (2*sxy)/3;
                    (4*syy)/9 - (2*sxx)/9 - (2*szz)/9, (4*sxx)/9 - (2*syy)/9 - (2*szz)/9, (4*szz)/9 - (2*syy)/9 - (2*sxx)/9, -(4*sxy)/3;
                    (2*sxy)/3, (2*sxy)/3, -(4*sxy)/3, (2*sxx)/3 + (2*syy)/3 - (4*szz)/3];
        end

        %------------------------------------------------------------------
        % von Mises stress
        function sVM = vonMisesStress(this,stress)
            J2 = this.stressInvariantJ2(stress);
            sVM = sqrt(3.0 * J2);
        end

        %------------------------------------------------------------------
        % von Mises stress gradient
        function dsVM = vonMisesStressGradient(this,stress)
            J2      = max(this.stressInvariantJ2(stress), 1e-12);
            dsVMdJ2 = 0.5*sqrt(3.0/J2);
            dJ2     = this.gradientJ2(stress);
            dsVM    = dsVMdJ2 * dJ2;
        end

        %------------------------------------------------------------------
        % von Mises stress hessian matrix
        function d2sVM = vonMisesStressHessian(this,stress)
            J2    = max(this.stressInvariantJ2(stress), 1e-12);
            dJ2   = this.gradientJ2(stress);
            d2J2  = this.hessianJ2();
            dsVMdJ2 = 0.5*sqrt(3.0/J2);
            d2sVM = dsVMdJ2 * d2J2 - (0.25 * sqrt(3.0/J2/J2/J2))*(dJ2 * dJ2');
        end

        %------------------------------------------------------------------
        % Get plane stress constitutive matrix
        function D = planeStressConstitutiveMatrix(~,D)
            % Get the sub-matrices
            D11 = [ D(1,1) , D(1,2) , D(1,4) ;
                    D(2,1) , D(2,2) , D(2,4) ;
                    D(4,1) , D(4,2) , D(4,4) ];
            D21 = [ D(3,1) ; D(3,2) ; D(3,4)];
            D12 = [ D(1,3) ; D(2,3) ; D(4,3)];
            D22 = D(3,3);
            % Compute the plane stress matrix
            Dt  = D11 - D12 * D21' / D22;
            % Assemble
            D = [ Dt(1,1) , Dt(1,2) , 0.0 , Dt(1,3) ;
                  Dt(2,1) , Dt(2,2) , 0.0 , Dt(2,3) ;
                   0.0    ,   0.0   , 1.0 ,   0.0   ;
                  Dt(3,1) , Dt(3,2) , 0.0 , Dt(3,3) ];
        end

        %% Strain utilities
        % Considering a 2D strain tensor: strain = [exx, eyy, ezz, 2exy]

        %------------------------------------------------------------------
        % Compute the elastic out-of-plane strain
        function elasticOutOfPlaneStrain(~,material,ip)
            if strcmp(ip.anm,'PlaneStress')
                nu = material.nu;
                ip.strain(3) = -nu/(1.0-nu)*(ip.strain(1)+ip.strain(2));
            end
        end

        %------------------------------------------------------------------
        % Strain tensor norm
        function normE = normTensor(~,strain)
            exx = strain(1);
            eyy = strain(2);
            ezz = strain(3);
            exy = strain(4) / 2.0;
            normE = exx * exx + eyy * eyy + ezz * ezz + 2.0 * exy * exy;
            normE = sqrt(normE);
        end

        %------------------------------------------------------------------
        % I1 strain invariant
        function I1 = strainInvariantI1(~,strain)
            I1 = strain(1) + strain(2) + strain(3);
        end

        %------------------------------------------------------------------
        % Volumetric strain
        function ev = volumetricStrain(this,strain)
            ev = this.strainInvariantI1(strain);
        end

        %------------------------------------------------------------------
        % Deviatoric strain
        function ed = deviatoricStrain(this,strain)
            ed = zeros(4,1);
            ev = this.volumetricStrain(strain);
            ed(1) = strain(1) - ev / 3.0;
            ed(2) = strain(2) - ev / 3.0;
            ed(3) = strain(3) - ev / 3.0;
            ed(4) = strain(4) / 2.0;
        end

        %------------------------------------------------------------------
        % Norm deviatoric strain
        function ned = normDeviatoricStrain(this,strain)
            ed = this.deviatoricStrain(strain);
            ned = ed(1)*ed(1) + ed(2)*ed(2) + ed(3)*ed(3) + 2*ed(4)*ed(4);
            ned = sqrt(ned);
        end

        % Gradient of the norm of the deviatoric strain
        function gradNormEd = gradientNormDeviatoricStrain(this,strain)
            ned = this.normDeviatoricStrain(strain);
            if ned > 0.0
                ed  = this.deviatoricStrain(strain);
                gradNormEd = ed / ned;
            else
                gradNormEd = zeros(4,1);
            end
        end

        % I2 strain invariant
        function I2 = strainInvariantI2(~,strain)
            exx = strain(1);
            eyy = strain(2);
            ezz = strain(3);
            exy = strain(4) / 2.0;
            I2 =  exx*eyy + eyy*ezz + exx*ezz - exy*exy;
        end

        %------------------------------------------------------------------
        % J2 strain invariant
        function J2 = strainInvariantJ2(this,strain)
            I1 = this.strainInvariantI1(strain);
            I2 = this.strainInvariantI2(strain);
            J2  = I1*I1/3.0 - I2;
            J2  = max(J2,0.0);
        end

        %------------------------------------------------------------------
        % Gradient of the J2 stress invariant
        function dI1 = gradientStrainInvariantI1(~)
            dI1 = [1.0;1.0;1.0;0.0];
        end

        %------------------------------------------------------------------
        % Gradient of the J2 strain invariant
        function dJ2 = gradientStrainInvariantJ2(~,strain)
            exx = strain(1);
            eyy = strain(2);
            ezz = strain(3);
            exy = strain(4) / 2.0;
            dJ2 = zeros(4,1);
            dJ2(1) = (2.0 * exx - eyy - ezz)/3.0;
            dJ2(2) = (2.0 * eyy - exx - ezz)/3.0;
            dJ2(3) = (2.0 * ezz - eyy - exx)/3.0;
            dJ2(4) = 2.0 * exy;
        end

        %------------------------------------------------------------------
        % Fourth order symmetric tensor
        function I4 = fourthOrderSymTensor(~)
            I4 = [ 1.0 , 0.0 , 0.0 , 0.0;
                   0.0 , 1.0 , 0.0 , 0.0;
                   0.0 , 0.0 , 1.0 , 0.0;
                   0.0 , 0.0 , 0.0 , 0.5 ];
        end
    end
    %% Static methods
    methods (Static)
        %------------------------------------------------------------------
        % Flag to determine if the material is elasto-plastic or not
        function flag = isElastoPlastic()
            flag = false;
        end
    end
end
