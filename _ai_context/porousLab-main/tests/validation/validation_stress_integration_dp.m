%% DESCRIPTION
%
% Validation script for mechanical constitutive models.
%
% Physics:
% * Mechanical (M)
%
% Authors:
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
%
close all; clear; clc;

%% MODEL DEFINITION

% Material properties
rock = PorousMedia('rock');
rock.mechanical    = 'elasticDruckerPrager';  % Mechanical constitutive law
rock.rho           = 2.0e+3;           % Density (kg/m3)
rock.Young         = 2.0e+7;           % Young modulus (Pa)
rock.nu            = 0.49;             % Poisson ratio
rock.cohesion      = 5.0e+4;           % Cohesion (Pa)
rock.frictionAngle = 20*pi/180;        % Friction angle (rad)

% Material parameters vector
mat = struct('porousMedia', rock);

% Integration point
ip = IntPoint([0.0,0.0], 1.0,  Material_M(mat));
ip.initializeMechanicalAnalysisModel('PlaneStrain');

%% RUN ANALYSIS

% Strain increment direction
dstrain0 = [0.0;  % exx
            0.0;  % eyy
            0.0;  % ezz
            1.0]; % gxy

mag = 1.0e-03;
ninc = 150;

idplot      = 4;
p  = zeros(ninc+1, 1);
q  = zeros(ninc+1, 1);

for i = 1:ninc
    ip.strain = ip.strainOld + mag * dstrain0;
    ip.stress = ip.mechanicalLaw();
    ip.stressOld = ip.stress;
    ip.strainOld = ip.strain;

    % Stress invariants
    p(i+1) = ip.constitutiveMdl.mechanical.hydrostaticStress(ip.stress);
    q(i+1) = ip.constitutiveMdl.mechanical.vonMisesStress(ip.stress);
end

%% POS-PROCESSING

%  Stress vs Total strain
figure;
grid on, box on, hold on;
plot(p, q, 'b-', 'LineWidth', 1.5);
xlabel('p');
ylabel('q');
set(gca, 'fontsize', 18, 'TickLabelInterpreter', 'latex');
