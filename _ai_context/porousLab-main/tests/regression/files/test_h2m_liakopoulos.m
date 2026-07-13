%% DESCRIPTION
%
% Reduced Liakopoulos problem based on
% examples/H2M/example_h2m_liakopoulos.m.
%
% Physics:
% * Two-phase flow hydro-mechanical (H2M)
%

%% MODEL

mdl = Model_H2M();
mdl.gravityOn = true;

%% MESH

Lx = 0.1;
Ly = 1.0;
Nx = 1;
Ny = 20;
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);
mdl.setMesh(node, elem);

%% MATERIALS

water = Fluid('water');
water.K = 2.2e9;

gas = IdealGas('gas');
gas.mu = 1.8e-5;
gas.M  = 0.028949;
gas.T  = 300;

rock = PorousMedia('rock');
rock.K                  = 4.5e-13;
rock.phi                = 0.2975;
rock.Ks                 = 1.0e+25;
rock.Slr                = 0.2;
rock.lambda             = 3.0;
rock.Young              = 1.3e+6;
rock.nu                 = 0.4;
rock.rho                = 2000.0;
rock.liqRelPermeability = 'Liakopoulos';
rock.gasRelPermeability = 'BrooksCorey';
rock.capillaryPressure  = 'Liakopoulos';
rock.setMinGasRelPermeability(1.0e-4);

mdl.setMaterial(rock, water, gas);

%% BOUNDARY AND INITIAL CONDITIONS

mdl.setDisplacementDirichletBCAtBorder('bottom', [NaN, 0.0]);
mdl.setDisplacementDirichletBCAtBorder('left',   [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('right',  [0.0, NaN]);

mdl.setPressureDirichletBCAtBorder('bottom', 101225.0);
mdl.setInitialPressureAtDomain(101225.0);

mdl.setGasPressureDirichletBCAtBorder('top',    101325.0);
mdl.setGasPressureDirichletBCAtBorder('bottom', 101325.0);
mdl.setInitialGasPressureAtDomain(101325.0);

mdl.preComputations();
for el = 1:mdl.nelem
    elem = mdl.element(el).type;
    for i = 1:elem.nIntPoints
        elem.intPoint(i).stressOld = [101325; 101325; 0.0; 0.0];
    end
end

%% PROCESS

ti        = 0.001;
dt        = 0.001;
tf        = 30.0;
dtmax     = 10.0;
dtmin     = 0.0001;
adaptStep = true;

anl = Anl_Transient("Newton");
anl.setUpTransientSolver(ti, dt, tf, dtmax, dtmin, adaptStep);
anl.setScaleLinearSystem(true);
anl.maxIter = 15;
anl.echo = false;
anl.run(mdl);
