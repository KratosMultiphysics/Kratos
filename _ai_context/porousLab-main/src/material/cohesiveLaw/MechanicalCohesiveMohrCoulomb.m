%% MechanicalCohesiveMohrCoulomb Class
% This class implements a Mohr-Coulomb cohesive law for mechanical
% discontinuities. It evaluates elastic trial stresses and applies the
% shear and tension cut-off return rules when the interface reaches the
% cohesive strength.
%
%% Methods
% * *eval*: Computes the interface stress vector and tangent matrix for the
%           given material and integration point.
% * *isElastoPlastic*: Static method that indicates that the material is
%                      elasto-plastic.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class Definition
classdef MechanicalCohesiveMohrCoulomb < MechanicalCohesiveLinearElastic  
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MechanicalCohesiveMohrCoulomb()
            this = this@MechanicalCohesiveLinearElastic();
            this.nstVar = this.nstVar + 1;  % Slip tendency
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        % Compute the stress vector and the constitutive matrix
        function [stress,Dt] = eval(this,material,ip)

            % Constitutive matrix
            Dt = this.elasticConstitutiveMatrix(material,ip);

            % Trial stress vector
            stress =  ip.stressOld;
            if norm(ip.strain) > 0.0
                stress = stress + Dt * (ip.strain - ip.strainOld);
            end
            ts = stress(1);
            tn = stress(2);

            % Get material parameters
            kn       = material.normalStiffness;
            ks       = material.shearStiffness;
            phi      = material.frictionAngle; tanPhi = tan(phi);
            psi      = material.dilationAngle; tanPsi = tan(psi);
            c        = material.cohesion;
            tncutoff = material.tensionCutOff;

            % Current shear strength
            tsr = c - tn * tanPhi;

            % Mohr-Coulomb yield criterion
            fMC = abs(ts) - tsr;

            % Evaluate the tension-cut-off surface
            fTC = tn - tncutoff;

            % Initialize the increment of the plastic strain vector
            dEp = zeros(ip.nVar,1);
           
            % Slip tendency parameter for an elastic case
            if tsr ~= 0.0
                ST = abs(ts) / tsr;
            else
                if abs(ts) > 1.0e-12
                    ST = 1.0;
                else
                    ST = 0.0;
                end
            end

            % Check if it is not an elastic step
            if (fMC > 1.0e-15) || (fTC > 1.0e-15)

                % Shear stress for which both surfaces intersects
                ts_P = sign(ts)*abs(c - tncutoff*tanPhi);

                % Denominator 
                keq = ks + kn * tanPhi * tanPsi;

                % Trial normal stress on the MC yield function
                lambdaMC = fMC / keq;
                tnMC = tn - kn * lambdaMC * tanPsi;

                if tnMC <= tncutoff % Return to the MC

                    % Vector normal to the MC plastic potential surface
                    n = [sign(ts); tanPsi];

                    % Plastic strain increment
                    dEp = lambdaMC * n;

                    % Traction vector
                    stress(1) = ts - ks * dEp(1);
                    stress(2) = tn - kn * dEp(2);

                    % Tangent constitutive matrix
                    Dt(1,1) = kn * ks * tanPhi * tanPsi / keq;
                    Dt(1,2) = - kn * ks * sign(ts) * tanPhi / keq;
                    Dt(2,1) = - kn * ks * sign(ts) * tanPsi / keq;
                    Dt(2,2) = kn * ks / keq;

                    % Slip tendency parameter
                    ST = 1.0;

                elseif (abs(ts) <= abs(ts_P)) % Return to the Cut-off

                    % Plastic multiplier
                    lambdaTC = fTC / kn;

                    % Vector normal to the MC plastic potential surface
                    n = [0.0; 1.0];

                    % Plastic strain increment
                    dEp = lambdaTC * n;

                    % Traction vector
                    stress = [ts; tncutoff];

                    % Tangent constitutive matrix
                    Dt = [ks  , 0.0;
                          0.0 , 0.0];

                    % Slip tendency parameter
                    if abs(ts_P) < 1.0e-12
                        ST = 1.0;
                    else
                        ST = abs(ts) / (c - tncutoff*tanPhi);
                    end

                else % Return to the intersection point

                    % Plastic multipliers
                    lambdaMC = (fMC - fTC*tanPhi)/ks;
                    lambdaTC = fTC / kn - lambdaMC * tanPsi;

                    % Vector normal to each plastic potential surface
                    nMC = [sign(ts); tanPsi];
                    nTC = [0.0; 1.0];

                    % Plastic strain increment
                    dEp = lambdaMC * nMC + lambdaTC * nTC;

                    % Traction vector
                    stress = [ts_P; tncutoff];

                    % Tangent constitutive matrix
                    Dt = zeros(ip.nVar,ip.nVar);

                    % Slip tendency parameter
                    ST = 1.0;
                end

            end

            % Update the current plastic strain at the integration point
            ip.plasticstrain = ip.plasticstrainOld + dEp;

            % Update the current slip tendency ratio
            ip.statevar(2) = max(ST,ip.statevarOld(2));


        end 
    end
    %% Public methods
    methods (Static)
        %------------------------------------------------------------------
        % Flag to return that the material is not elasto-plastic
        function flag = isElastoPlastic()
            flag = true;
        end
        
    end
end
