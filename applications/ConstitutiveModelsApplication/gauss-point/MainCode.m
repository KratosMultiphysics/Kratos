function [] = MainCode()
close all
global M;

M = 1.5;


MakeThisFile('GaussPoint.csv')

function [] = MakeThisFile(XFILE);

rawData = csvread(XFILE);

Stress = rawData(:,2:7);
Strain = rawData(:,8:13);
ps = rawData(:,14);
pt = rawData(:,15);
pc = rawData(:,16);

figure(1)
PlotYieldSurface( ps(1), pt(1), pc(1), 'k');
PlotYieldSurface( ps(end), pt(end), pc(end), 'g-.');

% Triaxial plane
p = -sum(Stress(:,1:3)')/3.0;
J = Stress(:,1)-Stress(:,2);
plot(p, J, 'b')
size(p)

figure(2)
AxialStrain= Strain(:,2);
plot( -AxialStrain, J);

figure(3)
VolStrain = -sum(Strain(:,1:3)');
semilogx( p, VolStrain);

function PlotYieldSurface( ps, pt, pc, SPEC)

global M;

pp = linspace(0, -pc, 400);
jj = M * sqrt(-pp.*(pp+pc));

plot( pp, M*pp, 'r-.')
hold on
plot(pp+pt, jj, SPEC, 'linewidth', 2.0)


pp = linspace(0, -ps, 400);
jj = M * sqrt(-pp.*(pp+ps));
plot(pp, jj, SPEC, 'linewidth', 1.0)