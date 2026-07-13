%% DESCRIPTION
%
% Block crossed by strong discontinuity.
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

mdl.equivalentContinuum = true;
mdl.condenseEnrDofs = true;

%% MESH

% Create mesh
Lx = 5.0;  % Horizontal dimension (m)
Ly = 3.0;  % Vertical dimension (m)
Nx = 5;    % Number of elements in the x-direction
Ny = 3;    % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);

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
mdl.setPressureDirichletBCAtBorder('left', 0.0);
mdl.setPressureDirichletBCAtBorder('right', 10.0);

%% DISCONTINUITIES

% Create discontinuities
Dx = [1.0; 4.0];  % X-coordinates of polyline defining the fracture
Dy = [1.1; 1.9];  % Y-coordinates of polyline defining the fracture
fracture = Discontinuity([Dx, Dy], true);

% Set fracture material properties
fracture.liquidFluid = water;
fracture.initialAperture = 1.0e-3;

% Add fractures to model
mdl.addPreExistingDiscontinuities(fracture);

%% PROCESS

% Run analysis
anl = Anl_Linear();
anl.run(mdl);

%% POST-PROCESS

% Print results to command window
mdl.printResults();

% Plot contours
mdl.plotField('Pressure');
hold on;
fracture.plotIntersectedGeometry();
