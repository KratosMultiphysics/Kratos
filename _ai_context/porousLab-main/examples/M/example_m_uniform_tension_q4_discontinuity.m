%% DESCRIPTION
%
% Uniform tension problem with internal interface.
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
% mdl.condenseEnrDofs   = false;
% mdl.subDivIntegration = true;
mdl.symmetricSDAEFEM  = false;
mdl.useNodalEnrDofs = true;

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
rock.Young = 10.0e+9;  % Young modulus (Pa)
rock.nu    = 0.3;      % Poisson ratio

% Set materials to model
mdl.setMaterial(rock);

%% BOUNDARY CONDITIONS

mdl.setDisplacementDirichletBCAtBorder('bottom', [NaN, 0.0]);
mdl.setDisplacementDirichletBCAtPoint([0.0, 0.0],[0.0, 0.0])
mdl.addLoadAtBorder('top', 2, 1.0e+3);

%% DISCONTINUITIES

% Create discontinuities 
Xd = [ 0.0 , 0.5*Ly;
       Lx  , 0.5*Ly];
fracture = Discontinuity(Xd, true);

% Set fracture material properties
fracture.cohesiveLaw     = 'elastic';
fracture.shearStiffness  = 1.0;       % Pa/m
fracture.normalStiffness = 1.0e12;    % Pa/m

% Add fractures to model
discontinuityData = struct('addTangentialStretchingMode', false, 'addNormalStretchingMode', false, 'addRelRotationMode', true);
mdl.addPreExistingDiscontinuities(fracture, discontinuityData);

%% PROCESS

% Run analysis
anl = Anl_Linear();
anl.run(mdl);

%% POST-PROCESS

% Plot model
mdl.plotField('Sy');
hold on;
fracture.plotIntersectedGeometry();

mdl.plotFieldAlongDiscontinuiy('Sn',1)
