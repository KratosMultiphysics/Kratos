%% DESCRIPTION
%
% Block crossed by strong discontinuity.
% Reference:
% Dias-da-Costa, D. (2013, June 21). 
% Finite elements with embedded discontinuities in the scope of the
% discrete crack approach [Lecture notes]. Thom Conference Room.
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
mdl.condenseEnrDofs = false;
mdl.useNodalEnrDofs = true;
mdl.symmetricSDAEFEM  = true;

%% MESH

% Create mesh
Lx = 2.0e-3;  % Horizontal dimension (m)
Ly = 2.0e-3;  % Vertical dimension (m)
Nx = 1;    % Number of elements in the x-direction
Ny = 1;    % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);

% Set mesh to model
mdl.setMesh(node,elem);

mdl.t = 1.0e-3;

%% MATERIALS

% Create porous media
rock = PorousMedia('rock');   
rock.Young = 1.0e+14;  % Young modulus (Pa)
rock.nu    = 0.0;      % Poisson ratio

% Set materials to model
mdl.setMaterial(rock);

%% BOUNDARY CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('bottom', [0.0, 0.0]);

% Loads
mdl.addLoadAtPoint([0.0,2.0], [-0.5, 1.5]);    
% mdl.addLoadAtPoint([2.0,2.0], [ 0.0, 2.0]);    
% mdl.addLoadAtPoint([0.0,2.0], [0.0, 1.0]);   

%% DISCONTINUITIES

% Create discontinuities 
Dx = [0.00; 2.00e-3];  % X-coordinates of polyline defining the fracture
Dy = [0.25e-3; 1.75e-3];  % Y-coordinates of polyline defining the fracture
% Dy = [1.0e-3; 1.0e-3];  % Y-coordinates of polyline defining the fracture
fracture = Discontinuity([Dx, Dy], true);

% Set fracture material properties
fracture.cohesiveLaw     = 'elastic';
fracture.initialAperture = 0.0;
fracture.shearStiffness  = 1.0e9;
fracture.normalStiffness = 1.0e9;

% Add fractures to model
discontinuityData = struct('addTangentialStretchingMode', false, ...
                           'addRelRotationMode', true);
mdl.addPreExistingDiscontinuities(fracture, discontinuityData);

%% PROCESS

% Run analysis
anl = Anl_Linear();
anl.run(mdl);

%% POST-PROCESS

% Print results to command window
mdl.printResults();

% Plot model
mdl.plotField('Model');
hold on;
fracture.plotIntersectedGeometry();
mdl.plotFieldAlongDiscontinuiy('Sn',1)
