%% DESCRIPTION
%
% Kueper and Frind problem using liquid pressure (Pl) and gas pressure (Pg)
% as primary variables.
%
% Physics:
% * Two-phase hydraulic (H2)
%
% Authors:
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
%
%% MODEL

% Create model
mdl = Model_H2();

% Set model options
mdl.gravityOn    = true;

%% MESH

% Create mesh
Lx = 0.7;  % Horizontal dimension (m)
Ly = 0.5;  % Vertical dimension (m)
Nx = 56;   % Number of elements in the x-direction
Ny = 40;   % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create fluids
water   = Fluid('water');
gas     = Fluid('gas');
gas.rho = 1.630e+3;  % Density (kg/m3)
gas.mu  = 0.900e-3;  % Viscosity (Pa*s)

% Create porous media
%                          | K(m2)   | phi | biot |  Ks   | Slr  | Sgr |  Pb  | lambda  | LiqRelPerm |  GasRelPerm  | capPressure
sand1 = PorousMedia('sand1', 5.04e-10, 0.40, 1.0,  1.0e+25, 0.078, 0.0, 369.73,  3.86,  'BrooksCorey', 'BrooksCorey','BrooksCorey');
sand2 = PorousMedia('sand2', 2.05e-10, 0.39, 1.0,  1.0e+25, 0.069, 0.0, 434.45,  3.51,  'BrooksCorey', 'BrooksCorey','BrooksCorey');
sand3 = PorousMedia('sand3', 5.26e-11, 0.39, 1.0,  1.0e+25, 0.098, 0.0, 1323.95, 2.49,  'BrooksCorey', 'BrooksCorey','BrooksCorey');
sand4 = PorousMedia('sand4', 8.19e-12, 0.41, 1.0,  1.0e+25, 0.189, 0.0, 3246.15, 3.30,  'BrooksCorey', 'BrooksCorey','BrooksCorey');

% Compute element centroids to set material IDs
nelem = Nx * Ny;
Xc = zeros(nelem, 2);
for i = 1:nelem
    xcoord = mdl.NODE(mdl.ELEM{i}, 1);
    ycoord = mdl.NODE(mdl.ELEM{i}, 2);
    xcentr = sum(xcoord) / 4;  % Considering linear quad elements
    ycentr = sum(ycoord) / 4;  % Considering linear quad elements
    Xc(i, :) = [xcentr, ycentr];
end

% Sand 1
mdl.matID = ones(nelem,1);

% Sand 2
reg = isInsideRectangle(Xc, [0.10,0.15], [0.60,0.20]); mdl.matID(reg==1) = 2;

% Sand 3
reg = isInsideRectangle(Xc, [0.10,0.20], [0.25,0.30]); mdl.matID(reg==1) = 3;
reg = isInsideRectangle(Xc, [0.35,0.20], [0.60,0.25]); mdl.matID(reg==1) = 3;
reg = isInsideRectangle(Xc, [0.20,0.35], [0.50,0.40]); mdl.matID(reg==1) = 3;

% Sand 4 
reg = isInsideRectangle(Xc, [0.00,0.00], [0.70,0.05]); mdl.matID(reg==1) = 4;
reg = isInsideRectangle(Xc, [0.05,0.05], [0.20,0.15]); mdl.matID(reg==1) = 4;
reg = isInsideRectangle(Xc, [0.50,0.05], [0.65,0.15]); mdl.matID(reg==1) = 4;
reg = isInsideRectangle(Xc, [0.05,0.15], [0.10,0.40]); mdl.matID(reg==1) = 4;
reg = isInsideRectangle(Xc, [0.20,0.10], [0.45,0.15]); mdl.matID(reg==1) = 4;
reg = isInsideRectangle(Xc, [0.60,0.15], [0.65,0.40]); mdl.matID(reg==1) = 4;
reg = isInsideRectangle(Xc, [0.35,0.25], [0.60,0.30]); mdl.matID(reg==1) = 4;

% Set materials to model
mdl.setMaterial([sand1, sand2, sand3, sand4], water, gas);

%% BOUNDARY AND INITIAL CONDITIONS

% Set prescribed pressures to top corners
mdl.setPressureDirichletBCAtPoint([0.0, Ly], 0.0);
mdl.setPressureDirichletBCAtPoint([Lx, Ly], 0.0);
mdl.setGasPressureDirichletBCAtPoint([0.0, Ly], 369.73);
mdl.setGasPressureDirichletBCAtPoint([Lx, Ly], 369.73);

% Set prescribed gas pressure to infiltration zone
tol = 1.0e-4;
reg = find(isInsideRectangle(mdl.NODE, [0.3-tol,0.5-tol], [0.4+tol,0.5+tol]));
for i = 1:length(reg)
    mdl.setPressureDirichletBCAtNode(reg(i), 215.063);
    mdl.setGasPressureDirichletBCAtNode(reg(i), 639.35);
end

% Set initial conditions
mdl.setInitialPressureAtDomain(0.0);
mdl.setInitialGasPressureAtDomain(369.73);

%% PROCESS

% Analysis parameters
ti        = 1.0;    % Initial time
dt        = 1.0;    % Time step
tf        = 184.0;  % Final time
dtmax     = 1.0;    % Maximum time step
dtmin     = 0.001;  % Minimum time step
adaptStep = true;   % Adaptive step size

% Run analysis
anl = Anl_Transient("Newton");
anl.setUpTransientSolver(ti, dt, tf, dtmax, dtmin, adaptStep);
anl.maxIter = 10;
anl.run(mdl);

%% POST-PROCESS

% Plot contours
mdl.plotField('Model');
mdl.plotField('GasSaturation', [0.0, 1.0]);
