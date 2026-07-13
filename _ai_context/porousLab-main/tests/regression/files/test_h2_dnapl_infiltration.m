%% DESCRIPTION
%
% Reduced DNAPL infiltration problem based on
% examples/H2/example_h2_dnapl_infiltration.m.
%
% Physics:
% * Two-phase hydraulic (H2)
%

%% MODEL

mdl = Model_H2();
mdl.gravityOn = true;

%% MESH

Lx = 0.90;
Ly = 0.65;
Nx = 30;
Ny = 20;
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);
mdl.setMesh(node, elem);

%% MATERIALS

water     = Fluid('water');
dnapl     = Fluid('dnapl');
dnapl.rho = 1.460e+3;
dnapl.mu  = 0.900e-3;

mat1 = PorousMedia('rock');
mat1.K                  = 6.64e-11;
mat1.phi                = 0.4;
mat1.biot               = 1.0;
mat1.Ks                 = 1.0e+25;
mat1.Slr                = 0.1;
mat1.Sgr                = 0.0;
mat1.Pb                 = 755.0;
mat1.lambda             = 2.7;
mat1.liqRelPermeability = 'BrooksCorey';
mat1.gasRelPermeability = 'BrooksCorey';
mat1.capillaryPressure  = 'BrooksCorey';

mat2 = PorousMedia('rock');
mat2.K                  = 6.64e-13;
mat2.phi                = 0.4;
mat2.biot               = 1.0;
mat2.Ks                 = 1.0e+25;
mat2.Slr                = 0.1;
mat2.Sgr                = 0.0;
mat2.Pb                 = 755.0;
mat2.lambda             = 2.7;
mat2.liqRelPermeability = 'BrooksCorey';
mat2.gasRelPermeability = 'BrooksCorey';
mat2.capillaryPressure  = 'BrooksCorey';

nelem = Nx * Ny;
Xc = zeros(nelem, 2);
for i = 1:nelem
    xcoord = mdl.NODE(mdl.ELEM{i}, 1);
    ycoord = mdl.NODE(mdl.ELEM{i}, 2);
    Xc(i, :) = [sum(xcoord) / 4, sum(ycoord) / 4];
end

mdl.matID = ones(nelem, 1);
reg = isInsideRectangle(Xc, [0.3,0.325], [0.45,0.4875]);
mdl.matID(reg==1) = 2;

mdl.setMaterial([mat1, mat2], water, dnapl);

%% BOUNDARY AND INITIAL CONDITIONS

pci = 755.0;

for i = 1:mdl.nnodes
    pg = 7886.5 - 9810.0 * mdl.NODE(i,2);
    pl = pg - pci;
    mdl.setInitialGasPressureAtNode(i, pg);
    mdl.setInitialPressureAtNode(i, pl);

    if abs(mdl.NODE(i,1)) < 1.0e-12 || abs(mdl.NODE(i,1) - Lx) < 1.0e-12
        mdl.setPressureDirichletBCAtNode(i, pl);
        mdl.setGasPressureDirichletBCAtNode(i, pg);
    end
end

qginj = 0.3 * 0.075 / dnapl.rho;
tol = 1.0e-4;
reg = find(isInsideRectangle(mdl.NODE, [0.3-tol,Ly-tol], [0.6+tol,Ly+tol]));
for i = 1:length(reg)
    mdl.setGasPressureNeumannBCAtNode(reg(i), qginj/length(reg));
end

%% PROCESS

ti        = 1.0;
dt        = 1.0;
tf        = 5.0;
dtmax     = 20.0;
dtmin     = 0.01;
adaptStep = true;

anl = Anl_Transient("Newton");
anl.setUpTransientSolver(ti, dt, tf, dtmax, dtmin, adaptStep);
anl.maxIter = 10;
anl.echo = false;
anl.run(mdl);
