%% DESCRIPTION
%
% Block crossed by strong discontinuity. Stretching mode.
%
% Reference DOI: 10.1016/j.cma.2011.05.008
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
mdl.subDivIntegration = true;
mdl.symmetricSDAEFEM  = false;

%% MESH

% Create mesh
Lx = 2.0e-1;  % Horizontal dimension (m)
Ly = 2.0e-1;  % Vertical dimension (m)
Nx = 1;       % Number of elements in the x-direction
Ny = 1;       % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);

% Set mesh to model
mdl.setMesh(node,elem);

%% MATERIALS

% Create porous media
rock = PorousMedia('rock');   
rock.Young = 30.0e+9;  % Young modulus (Pa)
rock.nu    = 0.0;      % Poisson ratio

% Set materials to model
mdl.setMaterial(rock);

%% BOUNDARY CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('left', [0.0, 0.0]);
mdl.setDisplacementDirichletBCAtPoint([2.0e-1,0.0],[0.0,-1.0e-4]);
mdl.setDisplacementDirichletBCAtPoint([2.0e-1,2.0e-1],[0.0,1.0e-4]);

%% DISCONTINUITIES

% Create discontinuities 
Dx = [1.0e-1; 1.0e-1];  % X-coordinates of polyline defining the fracture
Dy = [0.0e-1; 2.0e-1];  % Y-coordinates of polyline defining the fracture
fracture = Discontinuity([Dx, Dy], true);

% Set fracture material properties
fracture.cohesiveLaw     = 'elastic';
fracture.initialAperture = 0.0;
fracture.shearStiffness  = 0.0;    % Pa/m
fracture.normalStiffness = 0.0;    % Pa/m

% Add fractures to model
discontinuityData = struct('addTangentialStretchingMode', true, 'addRelRotationMode', true, 'addNormalStretchingMode', false);
mdl.addPreExistingDiscontinuities(fracture, discontinuityData);

%% PROCESS

% Run analysis
anl = Anl_Linear();
anl.run(mdl);

%% POST-PROCESS

% Print results to command window
mdl.printResults();

% Print stresses
fprintf("\n\n IP      Sxx          Syy          Szz          Sxy\n")
elem = mdl.element(1).type;
for i = 1:elem.nIntPoints
    stress = elem.intPoint(i).stress;
    fprintf("%2d   %+.4e  %+.4e  %+.4e  %+.4e \n",i,stress(1),stress(2),stress(3),stress(4));
end

% Plot model
mdl.plotField('Ux');
hold on;
fracture.plotIntersectedGeometry();
