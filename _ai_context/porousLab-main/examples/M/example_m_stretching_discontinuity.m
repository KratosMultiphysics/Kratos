%% DESCRIPTION
%
% Block crossed by strong discontinuity. Stretching mode.
%
% Reference DOI: 10.1016/j.cma.2009.07.013
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
mdl.subDivIntegration = true;
mdl.useNodalEnrDofs = true;

%% MESH

% Create mesh
Lx = 2.0e-3;  % Horizontal dimension (m)
Ly = 1.0e-3;  % Vertical dimension (m)
Nx = 1;       % Number of elements in the x-direction
Ny = 1;       % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);

% Model thickness
mdl.t = 1.0e-3;

% Set mesh to model
mdl.setMesh(node,elem);

%% MATERIALS

% Create porous media
rock = PorousMedia('rock');   
rock.Young = 30.0e+6;  % Young modulus (kPa)
rock.nu    = 0.0;      % Poisson ratio

% Set materials to model
mdl.setMaterial(rock);

%% BOUNDARY CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('bottom', [0.0, 0.0]);
mdl.setDisplacementDirichletBCAtBorder('top',    [NaN, 0.0]);

% Loads
mdl.addLoadAtPoint([0.0,1.0], [ 13.49,0.0]);      
mdl.addLoadAtPoint([2.0,1.0], [-13.49,0.0]);    

%% DISCONTINUITIES

% Create discontinuities 
Dx = [0.0; 2.0e-3];     % X-coordinates of polyline defining the fracture
Dy = [0.5e-3; 0.5e-3];  % Y-coordinates of polyline defining the fracture
fracture = Discontinuity([Dx, Dy], true);

% Set fracture material properties
fracture.cohesiveLaw     = 'elastic';
fracture.initialAperture = 0.0;
fracture.shearStiffness  = 1.0e10;    % Pa/m
fracture.normalStiffness = 1.0e10;    % Pa/m

% Add fractures to model
discontinuityData = struct('addTangentialStretchingMode', true, 'addNormalStretchingMode', false, 'addRelRotationMode', true);
mdl.addPreExistingDiscontinuities(fracture, discontinuityData);

%% PROCESS

% Run analysis
anl = Anl_Linear();
anl.run(mdl);

%% POST-PROCESS

% Print results to command window
mdl.printResults();

% Plot model
mdl.plotField('Ux');
hold on;
fracture.plotIntersectedGeometry();
