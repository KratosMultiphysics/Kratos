%% DESCRIPTION
%
% Reduced CO2 injection problem based on
% examples/H2/example_h2_co2_injection.m.
%
% Physics:
% * Two-phase hydraulic (H2)
%

%% MODEL

mdl = Model_H2();
mdl.isAxisSymmetric = true;
mdl.gravityOn       = true;

%% MESH

Lx = 200.0;
Ly = 6.0;
Nx = 8;
Ny = 12;
type = 'ISOQ4';
quadDistrX = true;
quadDistrY = false;
[node, elem] = regularMesh(Lx, Ly, Nx, Ny, [], [], type, quadDistrX, quadDistrY);
mdl.setMesh(node, elem);

%% MATERIALS

brine     = Fluid('brine');
brine.rho = 1173.0;
brine.mu  = 1.252e-3;

co2     = Fluid('co2');
co2.rho = 848.0;
co2.mu  = 8.100e-5;

aquifer = PorousMedia('aquifer');
aquifer.K                  = 3.0e-12;
aquifer.phi                = 0.26;
aquifer.biot               = 1.0;
aquifer.Ks                 = 1.0e+25;
aquifer.Slr                = 0.35;
aquifer.Sgr                = 0.0;
aquifer.Pb                 = 1.0e+4;
aquifer.lambda             = 2.0;
aquifer.liqRelPermeability = 'BrooksCorey';
aquifer.gasRelPermeability = 'BrooksCorey';
aquifer.capillaryPressure  = 'BrooksCorey';

mdl.setMaterial(aquifer, brine, co2);

%% BOUNDARY AND INITIAL CONDITIONS

depth1 = 1500.0;
depth0 = depth1 + Ly;
grav = 9.806;
pci = 1.0e4;

for i = 1:mdl.nnodes
    h = depth0 - mdl.NODE(i,2);
    pl = grav * brine.rho * h;
    pg = pci + pl;

    mdl.setInitialPressureAtNode(i, pl);
    mdl.setInitialGasPressureAtNode(i, pg);

    if abs(mdl.NODE(i,1) - Lx) < 1.0e-12
        mdl.setPressureDirichletBCAtNode(i, pl);
        mdl.setGasPressureDirichletBCAtNode(i, pg);
    end
end

day = 60 * 60 * 24;
qinj = (1.0 / day) * Ly / (Ny + 1);
mdl.setGasPressureNeumannBCAtBorder('left', qinj);

%% PROCESS

ti        = 0.005*day;
dt        = 0.005*day;
tf        = 0.05*day;
dtmax     = 0.01*day;
dtmin     = 1.0e-4*day;
adaptStep = true;

anl = Anl_Transient("Newton");
anl.setUpTransientSolver(ti, dt, tf, dtmax, dtmin, adaptStep);
anl.maxIter = 10;
anl.echo = false;
anl.run(mdl);
