%% DESCRIPTION
%
% Dense non-liquid phase infiltration problem using liquid pressure (Pl)
% and gas pressure (Pg) as primary variables.
%
% References:
% * Benisch et al (2013). The coupled OpenGeoSys-eclipse simulator for simulation of CO2 storage–code comparison for fluid flow and geomechanical processes. Energy Procedia, 37:3663-3671.
% * Graupner et al (2011). The coupled simulator ECLIPSE–OpenGeoSys for the simulation of CO2 storage in saline formations. Energy Procedia, 4:3794-3800.
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
mdl.isAxisSymmetric = true;
mdl.gravityOn       = true;

%% MESH

% Create mesh
Lx = 200.0;          % Horizontal dimension (m)
Ly = 6.0;            % Vertical dimension (m)
Nx = 30;             % Number of elements in the x-direction
Ny = 60;             % Number of elements in the y-direction
type = 'ISOQ4';      % Type of finite element
quadDistrX = true;   % Flag for quadratic distribution in x-direction 
quadDistrY = false;  % Flag for quadratic distribution in y-direction 
[node, elem] = regularMesh(Lx, Ly, Nx, Ny, [], [], type, quadDistrX, quadDistrY);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create fluids
brine     = Fluid('brine');
brine.rho = 1173.0;    % Density (kg/m3)
brine.mu  = 1.252e-3;  % Viscosity (Pa*s)

co2     = Fluid('co2');
co2.rho = 848.0;     % Density (kg/m3)
co2.mu  = 8.100e-5;  % Viscosity (Pa*s)

% Create porous media
aquifer = PorousMedia('aquifer');
aquifer.K                  = 3.0e-12;        % Intrinsic permeability (m2)
aquifer.phi                = 0.26;           % Porosity
aquifer.biot               = 1.0;            % Biot's coefficient
aquifer.Ks                 = 1.0e+25;        % Solid bulk modulus (Pa)
aquifer.Slr                = 0.35;           % Residual liquid saturation
aquifer.Sgr                = 0.0;            % Residual gas saturation
aquifer.Pb                 = 1.0e+4;         % Gas-entry pressure
aquifer.lambda             = 2.0;            % Curve-fitting parameter
aquifer.liqRelPermeability = 'BrooksCorey';  % Liquid relative permeability
aquifer.gasRelPermeability = 'BrooksCorey';  % Gas relative permeability
aquifer.capillaryPressure  = 'BrooksCorey';  % Saturation degree function

% Set materials to model
mdl.setMaterial(aquifer, brine, co2);

%% BOUNDARY AND INITIAL CONDITIONS

% Gas pressure (hydrostatic profile)
depth1 = 1500.0;       % Depth of the aquifer (top border)
depth0 = depth1 + Ly;  % Depth of the aquifer (bottom border)
grav = 9.806;          % Gravity acceleration (m/s2)
pci = 1.0e4;           % Initial capillary pressure (Pa)

for i = 1:mdl.nnodes
    h = depth0 - mdl.NODE(i,2);
    pl = grav * brine.rho * h;
    pg = pci + pl;

    % Set initial conditions
    mdl.setInitialPressureAtNode(i, pl);
    mdl.setInitialGasPressureAtNode(i, pg);

    % Fix pressures at right borders
    if ((abs(mdl.NODE(i,1) - Lx)) < 1.0e-12)
        mdl.setPressureDirichletBCAtNode(i, pl);
        mdl.setGasPressureDirichletBCAtNode(i, pg);
    end
end

% Set prescribed gas pressure to infiltration zone
day = 60 * 60 * 24;
qinj = (1.0 / day) * Ly / (Ny + 1);
mdl.setGasPressureNeumannBCAtBorder('left', qinj);

%% PROCESS

% Analysis parameters
ti        = 0.005*day;    % Initial time
dt        = 0.005*day;    % Time step
tf        = 50.00*day;    % Final time
dtmax     = 1.0*day;     % Maximum time step
dtmin     = 1.0e-3*day;  % Minimum time step
adaptStep = true;        % Adaptive step size

% Run analysis
anl = Anl_Transient("Newton");
anl.setUpTransientSolver(ti, dt, tf, dtmax, dtmin, adaptStep);
anl.maxIter = 10;
anl.run(mdl);

%% POST-PROCESS

% Plot contours
mdl.plotField('CapillaryPressure');
mdl.plotField('GasPressure');

% Plot graphs
Xi = [0.0, 0.0]; Xf = [Lx, 0.0];
mdl.plotFieldAlongSegment('LiquidPressure', Xi, Xf, 500, 'x');
