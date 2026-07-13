%% DESCRIPTION
%
% Physics:
% * Single-phase hydraulic (H)
%
% Authors:
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
%
%% MODEL

% Create model
mdl = Model_H();

%% MESH

% Create mesh
[node, elem] = regularMesh(2.0, 1.0, 2, 1);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create fluids
water = Fluid('water');
water.rho = 1.0e+3;  % Density (kg/m3)
water.mu  = 1.0e-3;  % Viscosity (Pa*s)
water.K   = 2.0e+9;  % Compressibility/Bulk modulus (1/Pa)

% Create porous media
rock = PorousMedia('rock');
rock.K    = 1.0194e-14;  % Intrinsic permeability (m2)
rock.phi  = 0.3;         % Porosity
rock.Ks   = 1.0e+12;     % Rock bulk modulus (Pa)
rock.biot = 0.6;         % Biot coefficient

% Set materials to model
mdl.setMaterial(rock, water);

%% BOUNDARY CONDITIONS

% Set Dirichlet boundary conditions
mdl.setPressureDirichletBCAtBorder('left', 0.0);
mdl.setPressureDirichletBCAtBorder('right', 10.0);

%% PROCESS

% Analysis parameters
ti        = 1.0;    % Initial time
dt        = 1.0;    % Time step
tf        = 100.0;  % Final time
dtmax     = 50.0;   % Maximum time step
dtmin     = 0.001;  % Minimum time step
adaptStep = true;   % Adaptive step size

% Run analysis
anl = Anl_Transient("Picard");
anl.setUpTransientSolver(ti, dt, tf, dtmax, dtmin, adaptStep);
anl.setRelativeConvergenceCriteria(true);
anl.echo = false;
anl.run(mdl);