%% DESCRIPTION
%
% Strip footing problem.
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
mdl.intOrder = 2;  % Integration quadrature order

%% MESH

% Load mesh
% load('MeshStripFooting');
% Create mesh
[node, elem] = regularMesh(5.0, 5.0, 30, 30, [], [], 'ISOQ4', 0.1, 0.92, 1.0, 0.7);
[node, elem] = convertToQuadraticMesh(node, elem);

% Set mesh to model
mdl.setMesh(node, elem);
mdl.resequenceNodes() ;

%% MATERIALS

% Create porous media
rock = PorousMedia('rock');
rock.mechanical    = 'druckerPrager';  % Constitutive law
rock.MCmatch       = 'planestrain';    % How Drucker-Prager surfaces matches Mohr-Coulomb
rock.Young         = 1.0e+7;           % Young modulus (kPa)
rock.nu            = 0.48;             % Poisson ratio
rock.cohesion      = 490.0;            % Cohesion (kPa)
rock.frictionAngle = 20*pi/180;        % Friction angle (rad)
rock.dilationAngle = 20*pi/180;        % Dilation angle (rad)
% rock.stressIntAlgorithm = 'alternative';

% Set materials to model
mdl.setMaterial(rock);

%% BOUNDARY CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('left',   [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('right',  [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('bottom', [NaN, 0.0]);

% Prandtl's limit pressure
Nq = exp(pi*tan(rock.frictionAngle)) * tan(pi/4 + rock.frictionAngle/2)^2;
Nc = (Nq - 1.0) * cot(rock.frictionAngle);
Plim = rock.cohesion * Nc;
mdl.addLoadAtBorder('top', 2, -Plim, [0.0 , 0.5]);

%% PROCESS

% Setup analysis
anl = Anl_NonlinearQuasiStatic('GeneralizedDisplacement');
anl.adjustStep    = true;
anl.increment     = 0.1;
anl.max_increment = 1.0;
anl.max_lratio    = 2.0;
anl.max_step      = 100;
anl.max_iter      = 20;
anl.trg_iter      = 4;
anl.tol           = 1.0e-4;

% Node and DOF used to plot Load Factor vs Displacement
ndId = mdl.closestNodeToPoint([0.0, 5.0]);
anl.setPlotDof(ndId, 2);

% Run analysis
anl.run(mdl);

%% POST-PROCESS

% Plot Load Factor vs Displacement
anl.plotCurves();

% Plot contours
mdl.plotField('PEMAG');
mdl.plotField('PEMAG',[0.0,0.01]);
mdl.plotField('Uy');
