%% DESCRIPTION
%
% Elastoplastic pressurized cylinder.
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
Lx = 1.0;  % Horizontal dimension (m)
Ly = 1.0;  % Vertical dimension (m)
Nx = 10;   % Number of elements in the x-direction
Ny = 10;   % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);

% Transform to cylindrical coordinates
ri = 0.1;  % Internal cylinder radius
re = 0.2;  % External cylinder radius
r = ri + node(:, 1) * (re - ri);
theta = node(:, 2) * (pi / 2);
node = [r .* cos(theta), r .* sin(theta)];

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create porous media
rock = PorousMedia('rock');
rock.mechanical = 'vonMises';  % Constitutive law
rock.Young      = 2.1e+11;     % Young modulus (Pa)
rock.nu         = 0.3;         % Poisson ratio
rock.sy0        = 2.4e+8;      % Initial yield stress (Pa)
rock.Kp         = 0.0;         % Plastic modulus (Pa)

% Set materials to model
mdl.setMaterial(rock);

%% BOUNDARY CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('left',   [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('bottom', [NaN, 0.0]);

% Radius of each node wrt to center at (0,0)
r  = sqrt((mdl.NODE(:,1)).^2 + (mdl.NODE(:,2)).^2);
sn = mdl.NODE(:,2) ./ r;
cs = mdl.NODE(:,1) ./ r;

% Loaded nodes
internalNodes = (r-ri) < eps;
nInternalNodes = sum(internalNodes);

% Internal pressure (Pa) and Force magnitude
pint = 192.0905814164710e+06;
F0 = (pint*mdl.t*ri*pi/2.0)/(nInternalNodes-1)/2.0;

% Occurrences of each node
allNodes = cell2mat(mdl.ELEM(:));
nodeCount = histcounts(allNodes, 1:(size(mdl.NODE,1)+1))';

% Apply internal pressure to internal face nodes
mdl.LOAD(internalNodes,:) = F0 * nodeCount(internalNodes) .* [cs(internalNodes) , sn(internalNodes)];

%% PROCESS

% Setup analysis
anl = Anl_NonlinearQuasiStatic('ArcLengthCylControl');
anl.adjustStep    = true;
anl.increment     = 0.1;
anl.max_increment = 0.1;
anl.max_lratio    = 2.0;
anl.max_step      = 100;
anl.max_iter      = 100;
anl.trg_iter      = 4;

% Node and DOF used to plot Load Factor vs Displacement
ndId = mdl.closestNodeToPoint([ri, 0.0]);
anl.setPlotDof(ndId, 1);

% Run analysis
anl.run(mdl);

%% POST-PROCESS

% Plot Load Factor vs Displacement
anl.plotCurves();

% Plot contours
mdl.plotField('S1');
mdl.plotField('Sr');
mdl.plotField('Sx');
