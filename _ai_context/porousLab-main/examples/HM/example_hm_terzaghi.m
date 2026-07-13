%% DESCRIPTION
%
% Terzaghi consolidation problem.
%
% Physics:
% * Single-phase flow hydro-mechanical (HM)
%
% Authors:
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
%
%% MODEL

% Create model
mdl = Model_HM();

%% MESH

% Create mesh
Lx = 0.1;  % Horizontal dimension (m)
Ly = 1.0;  % Vertical dimension (m)
Nx = 5;    % Number of elements in the x-direction
Ny = 50;   % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create fluids
water = Fluid('water');

% Create porous media
rock = PorousMedia('rock');
rock.K     = 1.15741e-12;   % Intrinsic permeability (m2)
rock.phi   = 0.3;           % Porosity
rock.Young = 1.0e+6;        % Young modulus (Pa)
rock.nu    = 0.3;           % Poisson ratio

% Set materials to model
mdl.setMaterial(rock, water);

%% BOUNDARY AND INITIAL CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('bottom', [0.0, 0.0]);
mdl.setDisplacementDirichletBCAtBorder('left',   [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('right',  [0.0, NaN]);

% Loads
mdl.addLoadAtBorder('top', 2, -1.0e4);

% Pressure
mdl.setPressureDirichletBCAtBorder('top', 0.0);

%% PROCESS

% Analysis parameters
ti = 1.0;    % Initial time
dt = 1.0;    % Time step
tf = 100.0;  % Final time

% Run analysis
anl = Anl_Transient("Newton");
anl.setUpTransientSolver(ti, dt, tf);
anl.run(mdl);

%% POST-PROCESS

% Plot contours
mdl.plotField('Pressure');

% Plot graphs
Xi = [0.0, 0.0]; Xf = [0.0, 1.0];
mdl.plotFieldAlongSegment('Pressure', Xi, Xf, 500, 'y');
mdl.plotFieldAlongSegment('Uy', Xi, Xf, 500, 'y');
