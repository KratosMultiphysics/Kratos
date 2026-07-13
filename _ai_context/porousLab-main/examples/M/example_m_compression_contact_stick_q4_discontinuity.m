%% DESCRIPTION
%
% Contact problem. Block with horizontal discontinuity in the stick
% condition.
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
mdl.isPlaneStress     = true;
mdl.condenseEnrDofs   = false;
% mdl.subDivIntegration = true;
mdl.symmetricSDAEFEM  = false;

%% MESH

% Create mesh
Lx = 1.0;      % Horizontal dimension (m)
Ly = 1.0;      % Vertical dimension (m)
Nx = 21;       % Number of elements in the x-direction
Ny = 21;       % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny, [], [], 'CST');

% Set mesh to model
mdl.setMesh(node,elem);

%% MATERIALS

% Create porous media
rock = PorousMedia('rock');   
rock.Young = 5.0e+9;  % Young modulus (Pa)
rock.nu    = 0.3;      % Poisson ratio

% Set materials to model
mdl.setMaterial(rock);

%% BOUNDARY CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('bottom', [0.0, 0.0]);
mdl.setDisplacementDirichletBCAtBorder('top', [0.0, -1.0e-3]);

%% DISCONTINUITIES

% Create discontinuities 
Xd = [ 0.0 , 0.51*Ly;
       Lx  , 0.51*Ly];
fracture = Discontinuity(Xd, true);

% Set fracture material properties
fracture.cohesiveLaw     = 'elastic';
fracture.shearStiffness  = 1.0;       % Pa/m
fracture.normalStiffness = 1.0e15;    % Pa/m

% Add fractures to model
discontinuityData = struct('addTangentialStretchingMode', false, 'addNormalStretchingMode', false, 'addRelRotationMode', false);
mdl.addPreExistingDiscontinuities(fracture, discontinuityData);

%% PROCESS

% Run analysis
anl = Anl_Linear();
anl.run(mdl);

%% POST-PROCESS

% Plot model
mdl.plotField('Uy');
hold on;
fracture.plotIntersectedGeometry();

mdl.plotFieldAlongDiscontinuiy('Sn',1)
