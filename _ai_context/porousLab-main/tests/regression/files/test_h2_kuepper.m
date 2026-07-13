%% DESCRIPTION
%
% Reduced Kueper and Frind problem based on
% examples/H2/example_h2_kuepper.m.
%
% Physics:
% * Two-phase hydraulic (H2)
%

%% MODEL

mdl = Model_H2();
mdl.gravityOn = true;

%% MESH

Lx = 0.7;
Ly = 0.5;
Nx = 14;
Ny = 10;
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);
mdl.setMesh(node, elem);

%% MATERIALS

water   = Fluid('water');
gas     = Fluid('gas');
gas.rho = 1.630e+3;
gas.mu  = 0.900e-3;

sand1 = PorousMedia('sand1', 5.04e-10, 0.40, 1.0, 1.0e+25, 0.078, 0.0, 369.73, 3.86, 'BrooksCorey', 'BrooksCorey', 'BrooksCorey');
sand2 = PorousMedia('sand2', 2.05e-10, 0.39, 1.0, 1.0e+25, 0.069, 0.0, 434.45, 3.51, 'BrooksCorey', 'BrooksCorey', 'BrooksCorey');
sand3 = PorousMedia('sand3', 5.26e-11, 0.39, 1.0, 1.0e+25, 0.098, 0.0, 1323.95, 2.49, 'BrooksCorey', 'BrooksCorey', 'BrooksCorey');
sand4 = PorousMedia('sand4', 8.19e-12, 0.41, 1.0, 1.0e+25, 0.189, 0.0, 3246.15, 3.30, 'BrooksCorey', 'BrooksCorey', 'BrooksCorey');

nelem = Nx * Ny;
Xc = zeros(nelem, 2);
for i = 1:nelem
    xcoord = mdl.NODE(mdl.ELEM{i}, 1);
    ycoord = mdl.NODE(mdl.ELEM{i}, 2);
    Xc(i, :) = [sum(xcoord) / 4, sum(ycoord) / 4];
end

mdl.matID = ones(nelem,1);
reg = isInsideRectangle(Xc, [0.10,0.15], [0.60,0.20]); mdl.matID(reg==1) = 2;
reg = isInsideRectangle(Xc, [0.10,0.20], [0.25,0.30]); mdl.matID(reg==1) = 3;
reg = isInsideRectangle(Xc, [0.35,0.20], [0.60,0.25]); mdl.matID(reg==1) = 3;
reg = isInsideRectangle(Xc, [0.20,0.35], [0.50,0.40]); mdl.matID(reg==1) = 3;
reg = isInsideRectangle(Xc, [0.00,0.00], [0.70,0.05]); mdl.matID(reg==1) = 4;
reg = isInsideRectangle(Xc, [0.05,0.05], [0.20,0.15]); mdl.matID(reg==1) = 4;
reg = isInsideRectangle(Xc, [0.50,0.05], [0.65,0.15]); mdl.matID(reg==1) = 4;
reg = isInsideRectangle(Xc, [0.05,0.15], [0.10,0.40]); mdl.matID(reg==1) = 4;
reg = isInsideRectangle(Xc, [0.20,0.10], [0.45,0.15]); mdl.matID(reg==1) = 4;
reg = isInsideRectangle(Xc, [0.60,0.15], [0.65,0.40]); mdl.matID(reg==1) = 4;
reg = isInsideRectangle(Xc, [0.35,0.25], [0.60,0.30]); mdl.matID(reg==1) = 4;

mdl.setMaterial([sand1, sand2, sand3, sand4], water, gas);

%% BOUNDARY AND INITIAL CONDITIONS

mdl.setPressureDirichletBCAtPoint([0.0, Ly], 0.0);
mdl.setPressureDirichletBCAtPoint([Lx, Ly], 0.0);
mdl.setGasPressureDirichletBCAtPoint([0.0, Ly], 369.73);
mdl.setGasPressureDirichletBCAtPoint([Lx, Ly], 369.73);

tol = 1.0e-4;
reg = find(isInsideRectangle(mdl.NODE, [0.3-tol,Ly-tol], [0.4+tol,Ly+tol]));
for i = 1:length(reg)
    mdl.setPressureDirichletBCAtNode(reg(i), 215.063);
    mdl.setGasPressureDirichletBCAtNode(reg(i), 639.35);
end

mdl.setInitialPressureAtDomain(0.0);
mdl.setInitialGasPressureAtDomain(369.73);

%% PROCESS

ti        = 1.0;
dt        = 1.0;
tf        = 10.0;
dtmax     = 1.0;
dtmin     = 0.001;
adaptStep = true;

anl = Anl_Transient("Newton");
anl.setUpTransientSolver(ti, dt, tf, dtmax, dtmin, adaptStep);
anl.maxIter = 10;
anl.echo = false;
anl.run(mdl);
