%% DESCRIPTION
%
% Reduced Buckley-Leverett problem based on
% examples/H2/example_h2_buckley_leverett.m.
%
% Physics:
% * Two-phase hydraulic (H2)
%

%% MODEL

mdl = Model_H2();

%% MESH

Lx = 1000.0;
Ly = 100.0;
Nx = 100;
Ny = 1;
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);
mdl.setMesh(node, elem);

%% MATERIALS

water = Fluid('water');
gas   = Fluid('gas');

rock = PorousMedia('rock');
rock.K                  = 1.0e-7;
rock.phi                = 0.2;
rock.Slr                = 0.2;
rock.Sgr                = 0.2;
rock.Pb                 = 0.0;
rock.lambda             = 2.0;
rock.liqRelPermeability = 'BrooksCorey';
rock.gasRelPermeability = 'BrooksCorey';
rock.capillaryPressure  = 'UMAT';

SlPcUMAT = [1.0, 0.2;
            0.0, 0.8];
rock.setUMATCapillaryPressureCurve(SlPcUMAT);

mdl.setMaterial(rock, water, gas);

%% BOUNDARY AND INITIAL CONDITIONS

mdl.setPressureDirichletBCAtBorder('left', 0.0);
mdl.setGasPressureDirichletBCAtBorder('left', 0.0);

mdl.setInitialPressureAtDomain(0.0);
mdl.setInitialGasPressureAtDomain(1.0);

%% PROCESS

day       = 60*60*24;
ti        = 0.001*day;
dt        = 0.1*day;
tf        = 5.0*day;
dtmax     = 1.0*day;
dtmin     = 0.1*day;
adaptStep = true;

anl = Anl_Transient("Newton");
anl.setUpTransientSolver(ti, dt, tf, dtmax, dtmin, adaptStep);
anl.echo = false;
anl.run(mdl);
