%% DESCRIPTION
%
% Slope stability problem.
%
% Physics:
% * Mechanical (M)
%
% Authors:
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
%
%% MODEL

% Create model
mdl = Model_M();

% Set model options
mdl.gravityOn = true;
mdl.intOrder  = 2; % Integration quadrature order

%% MESH

% Load mesh
load('MeshSlopeStabilityTransfiniteQuadratic');

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create porous media
rock = PorousMedia('rock');
rock.mechanical    = 'druckerPrager';  % Constitutive law
rock.MCmatch       = "planestrain";    % How Drucker-Prager surfaces matches Mohr-Coulomb
rock.rho           = 2.0;              % Density (g/cm3)
rock.Young         = 2.0e+4;           % Young modulus (kPa)
rock.nu            = 0.49;             % Poisson ratio
rock.cohesion      = 50.0;             % Cohesion (kPa)
rock.frictionAngle = 20.0*pi/180;      % Friction angle (rad)
rock.dilationAngle = 20.0*pi/180;      % Dilation angle

% Set materials to model
mdl.setMaterial(rock);

%% BOUNDARY CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('left',   [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('right',  [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('bottom', [0.0, 0.0]);

%% PROCESS

% Setup analysis
anl = Anl_NonlinearQuasiStatic('ArcLengthCylControl');
anl.adjustStep    = true;
anl.increment     = 0.1;
anl.max_increment = 0.5;
anl.max_lratio    = 10.0;
anl.max_step      = 10;
anl.max_iter      = 100;
anl.trg_iter      = 9;

% Node and DOF used to plot Load Factor vs Displacement
ndId = mdl.closestNodeToPoint([35.0, 40.0]);
anl.setPlotDof(ndId, 2);

% Run analysis
anl.echo = false;
anl.run(mdl);
