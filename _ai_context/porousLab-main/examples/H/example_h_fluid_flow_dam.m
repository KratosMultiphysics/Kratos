%% DESCRIPTION
%
% Fluid-flow through the foundation of a gravity dam.
%
% References:
% * Segura and Carol (2004).
% On zero-thickness interface elements for diffusion problems.
% Int J Numer Anal Methods Geomech, 28(9):947-962.
%
% Physics:
% * Single-phase hydraulic (H)
%
% Authors:
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
%
%% MODEL

% Create model
mdl = Model_H();

mdl.intOrder = 3;

%% MESH

% Create mesh
Lx = 24.0;  % Horizontal dimension (m)
Ly = 6.0;   % Vertical dimension (m)
Nx = 48;  % Number of elements in the x-direction
Ny = 12;  % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create fluids
water = Fluid('water');

% Create porous media
rock = PorousMedia('rock');
rock.K   = 1.0194e-14;  % Intrinsic permeability (m2)
rock.phi = 0.3;         % Porosity

% Set materials to model
mdl.setMaterial(rock, water);

%% BOUNDARY CONDITIONS

% Set Dirichlet boundary conditions
mdl.setPressureDirichletBCAtBorder('top', 120.0e3, [0.0, 8.0]);
mdl.setPressureDirichletBCAtBorder('top', 60.0e3,  [12.0, 24.0]);

%% DISCONTINUITIES

% Generate fractures
% FractureDataDamFoundation = generateRandomFractures(Lx, Ly, 20, 2*Lx/Nx, 'Seed', 12345);
FractureDataDamFoundation = generateParallelFractures(Lx, Ly, -pi/4.0, 1.0);
% FractureDataDamFoundation2 = generateParallelFractures(Lx, Ly,  pi/3.0, 3.0);
% FractureDataDamFoundation = [FractureDataDamFoundation1; FractureDataDamFoundation2];

% Number of discontinuities
nd = length(FractureDataDamFoundation);

% Create discontinuities
fractures(1,nd) = Discontinuity();
for i = 1:nd
    fractures(i) = Discontinuity(FractureDataDamFoundation{i}, true);
end

% Set fracture material properties
rng(123,'twister');
for i = 1:nd
    fractures(i).liquidFluid = water;
    % fractures(i).initialAperture = 1e-5 + (5e-4 - 1e-5) .* rand;
    fractures(i).initialAperture = 1e-4;
end

% Add fractures to model
mdl.addPreExistingDiscontinuities(fractures);

%% PROCESS

% Run analysis
anl = Anl_Linear();
anl.run(mdl);

%% POST-PROCESS

% Print results to command window
mdl.printResults();

% Plot model
mdl.plotField('Model');
colorbar off; hold on;
for i = 1:length(mdl.discontinuitySet)
    mdl.discontinuitySet(i).plotIntersectedGeometry();
end

% Plot contours
mdl.plotField('Pressure',[50e3,120e3]);
hold on;
for i = 1:length(mdl.discontinuitySet)
    mdl.discontinuitySet(i).plotIntersectedGeometry();
end

% Plot graphs
Xi = [0.0, Ly/2.0]; Xf = [Lx, Ly/2.0];
mdl.plotFieldAlongSegment('Pressure', Xi, Xf, 500, 'x');
