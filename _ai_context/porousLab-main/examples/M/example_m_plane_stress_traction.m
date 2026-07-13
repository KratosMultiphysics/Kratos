%% DESCRIPTION
%
% Uniform traction on a plate with isotropic damage model.
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

%% MESH

% Create mesh
Lx = 0.11;  % Horizontal dimension (m)
Ly = 0.04;  % Vertical dimension (m)
Nx = 22;    % Number of elements in the x-direction
Ny = 8;     % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create porous media
rock = PorousMedia('rock');
rock.mechanical          = 'isoDamage';  % Constitutive law
rock.Young               = 2.0e+10;      % Young modulus (Pa)
rock.nu                  = 0.0;          % Poisson ratio
rock.kappa               = 10.0;         % Ratio of tensile and compressive strength
rock.DamageThreshold     = 1.0e-4;       % Damage threshold
rock.FractureEnergyMode1 = 50.0;         % Fracture energy associated with mode 1 (N/m)

% Set materials to model
mdl.setMaterial(rock);

%% BOUNDARY CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('left', [0.0, 0.0]);

% Loads
mdl.addLoadAtBorder('right', 1, 2.0e+6);

%% PROCESS

% Setup analysis
anl = Anl_NonlinearQuasiStatic('ArcLengthCylControl');
anl.adjustStep = true;
anl.increment  = 0.01;
anl.max_lratio = 200.0;
anl.max_step   = 15;
anl.max_iter   = 100;
anl.trg_iter   = 4;

% Node and DOF used to plot Load Factor vs Displacement
ndId = mdl.closestNodeToPoint([Lx, 0.0]);
anl.setPlotDof(ndId, 1);

% Run analysis
anl.run(mdl);

%% POST-PROCESS

% Plot Load Factor vs Displacement
anl.plotCurves();

% Plot contours
mdl.plotField('Ux');
mdl.plotField('Sx');
mdl.plotField('Sy');
mdl.plotField('Sxy');
