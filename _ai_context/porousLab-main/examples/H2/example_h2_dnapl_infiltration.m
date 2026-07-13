%% DESCRIPTION
%
% Dense non-liquid phase infiltration problem using the two-phase flow
% formulation.
%
% References:
% * Wang et al (2015). A parallel finite element method for two-phase flow processes in porous media: OpenGeoSys with PETSc. Environ Earth Sci, 73:2269–2285.
% * Bastian et al (2007). Numerical simulation and experimental studies of unsaturated water flow in heterogeneous systems. Reactive Flows, Diffusion and Transport. https://doi.org/10.1007/978-3-540-28396-6_22
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
Lx = 0.90;  % Horizontal dimension (m)
Ly = 0.65;  % Vertical dimension (m)
Nx = 60;    % Number of elements in the x-direction
Ny = 40;    % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create fluids
water     = Fluid('water');
dnapl     = Fluid('dnapl');
dnapl.rho = 1.460e+3;  % Density (kg/m3)
dnapl.mu  = 0.900e-3;  % Viscosity (Pa*s)

% Create porous media
mat1 = PorousMedia('rock');
mat1.K                  = 6.64e-11;       % Intrinsic permeability (m2)
mat1.phi                = 0.4;            % Porosity
mat1.biot               = 1.0;            % Biot's coefficient
mat1.Ks                 = 1.0e+25;        % Solid bulk modulus (Pa)
mat1.Slr                = 0.1;            % Residual liquid saturation
mat1.Sgr                = 0.0;            % Residual gas saturation
mat1.Pb                 = 755.0;          % Gas-entry pressure
mat1.lambda             = 2.7;            % Curve-fitting parameter
mat1.liqRelPermeability = 'BrooksCorey';  % Liquid relative permeability
mat1.gasRelPermeability = 'BrooksCorey';  % Gas relative permeability
mat1.capillaryPressure  = 'BrooksCorey';  % Saturation degree function

mat2 = PorousMedia('rock');
mat2.K                  = 6.64e-13;       % Intrinsic permeability (m2)
mat2.phi                = 0.4;            % Porosity
mat2.biot               = 1.0;            % Biot's coefficient
mat2.Ks                 = 1.0e+25;        % Solid bulk modulus (Pa)
mat2.Slr                = 0.1;            % Residual liquid saturation
mat2.Sgr                = 0.0;            % Residual gas saturation
mat2.Pb                 = 755.0;          % Gas-entry pressure
mat2.lambda             = 2.7;            % Curve-fitting parameter
mat2.liqRelPermeability = 'BrooksCorey';  % Liquid relative permeability
mat2.gasRelPermeability = 'BrooksCorey';  % Gas relative permeability
mat2.capillaryPressure  = 'BrooksCorey';  % Saturation degree function

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

% Porous media material 1
mdl.matID = ones(nelem, 1);

% Porous media material 2
reg = isInsideRectangle(Xc, [0.3,0.325], [0.45,0.4875]);
mdl.matID(reg==1) = 2;

% Set materials to model
mdl.setMaterial([mat1, mat2], water, dnapl);

%% BOUNDARY AND INITIAL CONDITIONS

% Initial capillary pressure (Pa)
pci = 755.0; 

% Hydrostatic profile
for i = 1:mdl.nnodes
    pg = 7886.5 - 9810.0 * mdl.NODE(i,2);
    pl = pg - pci;
    mdl.setInitialGasPressureAtNode(i, pg);
    mdl.setInitialPressureAtNode(i, pl);

    % Fix pressure at lateral borders
    if ((abs(mdl.NODE(i,1)) < 1.0e-12) || ((abs(mdl.NODE(i,1) - Lx)) < 1.0e-12))
        mdl.setPressureDirichletBCAtNode(i, pl);
        mdl.setGasPressureDirichletBCAtNode(i, pg);
    end
end

% Set prescribed gas pressure to infiltration zone
qginj = 0.3 * 0.075 / dnapl.rho;
tol = 1.0e-4;
reg = find(isInsideRectangle(mdl.NODE, [0.3-tol,Ly-tol], [0.6+tol,Ly+tol]));
for i = 1:length(reg)
    mdl.setGasPressureNeumannBCAtNode(reg(i), qginj/length(reg));
end

%% PROCESS

% Analysis parameters
ti        = 1.0;    % Initial time
dt        = 1.0;    % Time step
tf        = 2000.0;  % Final time
dtmax     = 20.0;   % Maximum time step
dtmin     = 0.01;   % Minimum time step
adaptStep = true;   % Adaptive step size

% Run analysis
anl = Anl_Transient("Newton");
anl.setUpTransientSolver(ti, dt, tf, dtmax, dtmin, adaptStep);
anl.maxIter = 10;
anl.run(mdl);

%% POST-PROCESS

% Plot contours
mdl.plotField('LiquidSaturation');
mdl.plotField('GasSaturation');

% Plot graphs
Xi = [0.0, 0.0]; Xf = [Lx, Ly];
mdl.plotFieldAlongSegment('LiquidPressure', Xi, Xf, 500, 'x');
