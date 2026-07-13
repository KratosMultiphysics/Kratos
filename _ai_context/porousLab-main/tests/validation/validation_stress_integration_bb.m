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

wmax = 1.0E-3;

mat = struct('fluid',[],...
             'cohesiveLaw','elastic', ...
             'initialAperture',0.001, ...
             'normalStiffness',10000.0, ...
             'shearStiffness',10000.0,...
             'contactPenalization','bartonBandis',...
             'maximumClosure',wmax,...
             'frictionAngle',[],...
             'dilationAngle', [],...
             'cohesion', [],...
             'tensionCutOff',[],...
             'leakoff',[]);

% Integration point
ip = IntPoint([0.0,0.0], 1.0,  MaterialDiscontinuity_M(mat));
ip.initializeMechanicalAnalysisModel('Interface');

%% RUN ANALYSIS

% Normal strain
strain_n = linspace(wmax,-2*wmax,100000);

idplot      = 2;
stressplot  = zeros(length(strain_n), 1);
strainplot  = zeros(length(strain_n), 1);

for i = 1:length(strain_n)
    % ip.strain = ip.strainOld + mag * dstrain0;
    ip.strain = [0.0; strain_n(i)];
    ip.stress = ip.mechanicalLaw();
    ip.stressOld = ip.stress;
    ip.strainOld = ip.strain;

    % Store
    stressplot(i+1)  = ip.stress(idplot);
    strainplot(i+1)  = ip.strain(idplot);
end

%% POS-PROCESSING

%  Stress vs Total strain
figure;
grid on, box on, hold on;
plot(strainplot, stressplot, 'b-', 'LineWidth', 1.5);
xlabel('\delta_{n}');
ylabel('t_{n}');
set(gca, 'fontsize', 18, 'TickLabelInterpreter', 'latex');
