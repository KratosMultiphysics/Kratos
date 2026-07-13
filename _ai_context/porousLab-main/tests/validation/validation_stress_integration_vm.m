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
rock.mechanical = 'vonMises'; % Elastoplastic with von Mises criteria 
rock.Young      = 2.0e+8;     % Young modulus (Pa)
rock.nu         = 0.3;        % Poisson ratio
rock.sy0        = 1.0e+5;     % Initial yield stress (Pa)
rock.Kp         = 0.0;        % Plastic modulus (Pa)

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

mag = 1.0e-5;
ninc = 150;

idplot      = 4;
stressplot  = zeros(ninc+1, 1);
strainplot  = zeros(ninc+1, 1);
pstrainplot = zeros(ninc+1, 1);

for i = 1:ninc
    ip.strain = ip.strainOld + mag * dstrain0;
    ip.stress = ip.mechanicalLaw();
    ip.stressOld = ip.stress;
    ip.strainOld = ip.strain;

    % Store
    stressplot(i+1)  = ip.stress(idplot);
    strainplot(i+1)  = ip.strain(idplot);
    pstrainplot(i+1) = ip.plasticstrain(idplot);
end

%% POS-PROCESSING

%  Stress vs Total strain
figure;
grid on, box on, hold on;
tauy = rock.sy0 / sqrt(3.0);
plot(strainplot, stressplot/tauy, 'b-', 'LineWidth', 1.5);
xlabel('\gamma_{xy}');
ylabel('\tau_{xy}/\tau_Y');
set(gca, 'fontsize', 18, 'TickLabelInterpreter', 'latex');

%  Stress vs Plastic strain
figure;
grid on, box on, hold on;
plot(pstrainplot, stressplot/tauy, 'r-', 'LineWidth', 1.5);
xlabel('\gamma_{xy}^p');
ylabel('\tau_{xy}/\tau_Y');
set(gca, 'fontsize', 18, 'TickLabelInterpreter', 'latex');
