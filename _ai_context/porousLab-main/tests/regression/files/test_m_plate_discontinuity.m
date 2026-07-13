%% DESCRIPTION
%
% Block crossed by strong discontinuity.
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
mdl.isPlaneStress   = true;
mdl.condenseEnrDofs = true;

%% MESH

% Create mesh
[node, elem] = regularMesh(2.0, 2.0, 1, 1);

% Set mesh to model
mdl.setMesh(node,elem);

%% MATERIALS

% Create porous media
rock = PorousMedia('rock');   
rock.Young = 1.0e+8;  % Young modulus (kPa)
rock.nu    = 0.0;     % Poisson ratio

% Set materials to model
mdl.setMaterial(rock);

%% BOUNDARY CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('bottom', [0.0, 0.0]);

% Loads
mdl.addLoadAtPoint([0.0,2.0], [-0.5,1.5]);                            

%% DISCONTINUITIES

% Create discontinuities 
Dx = [0.00; 2.00];  % X-coordinates of polyline defining the fracture
Dy = [0.25; 1.75];  % Y-coordinates of polyline defining the fracture
fracture = Discontinuity([Dx, Dy], true);

% Set fracture material properties
fracture.cohesiveLaw     = 'elastic';
fracture.initialAperture = 0.0;
fracture.shearStiffness  = 1.0;
fracture.normalStiffness = 1.0;

% Add fractures to model
discontinuityData = struct('addStretchingMode', false, 'addRelRotationMode', true);
mdl.addPreExistingDiscontinuities(fracture, discontinuityData);

%% PROCESS

% Run analysis
anl = Anl_Linear();
anl.run(mdl);