%% DESCRIPTION
%
% Uniform traction on a plate with linear elastic model and quadratic finite element mesh.
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
[node, elem] = regularMesh(0.11, 0.04, 22, 8);
[node, elem] = convertToQuadraticMesh(node, elem);

% Set mesh to model
mdl.setMesh(node, elem);
mdl.resequenceNodes();

%% MATERIALS

% Create porous media
rock = PorousMedia('rock');
rock.mechanical = 'elastic';  % Constitutive law
rock.Young      = 2.0e+10;    % Young modulus (Pa)
rock.nu         = 0.0;        % Poisson ratio

% Set materials to model
mdl.setMaterial(rock);

%% BOUNDARY CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('left', [0.0, 0.0]);

% Loads
mdl.addLoadAtBorder('right', 1, 2.0e+6);

%% PROCESS

% Run analysis
anl = Anl_Linear();
anl.run(mdl);