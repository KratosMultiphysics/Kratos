%% DESCRIPTION
%
% Reduced McWhorter and Sunada problem based on
% examples/H2/example_h2_mcwhorter.m.
%
% Physics:
% * Two-phase hydraulic (H2)
%

%% MODEL

mdl = Model_H2();

%% MESH

[node, elem] = regularMesh(2.6, 0.5, 100, 1);
mdl.setMesh(node, elem);

%% MATERIALS

water = Fluid('water');
gas   = Fluid('gas');
gas.mu = 0.005;

rock = PorousMedia('rock');
rock.K                  = 1.0e-10;
rock.phi                = 0.15;
rock.Slr                = 0.02;
rock.Sgr                = 0.001;
rock.Pb                 = 5.0e+3;
rock.lambda             = 3.0;
rock.liqRelPermeability = 'BrooksCorey';
rock.gasRelPermeability = 'BrooksCorey';
rock.capillaryPressure  = 'BrooksCorey';

mdl.setMaterial(rock, water, gas);

%% BOUNDARY AND INITIAL CONDITIONS

pc = @(Sl) (rock.Pb * ((Sl-rock.Slr)/(1.0-rock.Slr-rock.Sgr))^(-1/rock.lambda));
pc_ini = pc(0.05);
pc_left = pc(0.8);
pg_ini = 100000.0;

mdl.setPressureDirichletBCAtBorder('left', pg_ini - pc_left);
mdl.setGasPressureDirichletBCAtBorder('left', pg_ini);

mdl.setInitialPressureAtDomain(pg_ini - pc_ini);
mdl.setInitialGasPressureAtDomain(pg_ini);

%% PROCESS

ti        = 0.001;
dt        = 0.001;
tf        = 20.0;
dtmax     = 5.0;
dtmin     = 0.0000001;
adaptStep = true;

anl = Anl_Transient("Newton");
anl.setUpTransientSolver(ti, dt, tf, dtmax, dtmin, adaptStep);
anl.maxIter = 20;
anl.echo = false;
anl.run(mdl);
