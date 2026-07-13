%% DESCRIPTION
%
% Block crossed by strong discontinuity acting as a barrier.
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
[node, elem] = regularMesh(10.0, 10.0, 20, 20);

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
rock.K   = 9.8e-16;  % Intrinsic permeability (m2)
rock.phi = 0.25;     % Porosity

% Set materials to model
mdl.setMaterial(rock, water);

%% BOUNDARY CONDITIONS

% Set Dirichlet boundary conditions
mdl.setPressureDirichletBCAtBorder('bottom', 0.0);
mdl.setPressureDirichletBCAtBorder('top', 10.0);

%% DISCONTINUITY

% Create discontinuities
Dx = [0.0; 10.0];  % X-coordinates of polyline defining the fracture
Dy = [8.0; 3.0];  % Y-coordinates of polyline defining the fracture
fracture = Discontinuity([Dx, Dy], true);

% Set fracture material properties
fracture.liquidFluid = water;
fracture.initialAperture = 1.0e-3;
fracture.leakoff = 1.0e-19;

% Add fractures to model
mdl.addPreExistingDiscontinuities(fracture);

%% PROCESS

% Analysis parameters
ti        = 1.0;     % Initial time
dt        = 1.0;     % Time step
tf        = 8000.0;  % Final time
dtmax     = 500.0;   % Maximum time step
dtmin     = 0.001;   % Minimum time step
adaptStep = true;    % Adaptive step size

% Initialize
anl = Anl_Transient("Newton");
anl.setUpTransientSolver(ti, dt, tf, dtmax, dtmin, adaptStep);
anl.run(mdl);

%% POST-PROCESS

% Print results to command window
mdl.printResults();

% Plot contours
mdl.plotField('Pressure');
hold on;
fracture.plotIntersectedGeometry();
