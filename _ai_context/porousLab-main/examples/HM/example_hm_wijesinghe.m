%% DESCRIPTION
%
% Wijesinghe problem.
%
% Physics:
% * Single-phase flow hydro-mechanical (HM)
%
% Authors:
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
%
%% MODEL

% Create model
mdl = Model_HM();

%% MESH

% Create mesh
[node, elem] = regularMesh(25.0, 1.0, 50, 21, [], [], 'ISOQ4', 0.0, 0.8, 0.5, 0.95);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create fluids
water = Fluid('water');

% Create porous media
rock = PorousMedia('rock');
rock.K     = 1.0e-21;        % Intrinsic permeability (m2)
rock.phi   = 0.001;         % Porosity
rock.Young = 60.0e+9;       % Young modulus (Pa)
rock.nu    = 0.0;           % Poisson ratio

% Set materials to model
mdl.setMaterial(rock, water);

%% BOUNDARY AND INITIAL CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('bottom', [0.0, 0.0]);
mdl.setDisplacementDirichletBCAtBorder('left',   [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('right',  [0.0, NaN]);

% Loads
mdl.addLoadAtBorder('top', 2, -50.0e6);

% Pressure
mdl.setPressureDirichletBCAtDomain(11.0e6);

%% DISCONTINUITY

% Polyline that defines the discontinuity
Xd = [ 0.0  , 0.5;
       25.0 , 0.5];
fracture = Discontinuity(Xd, true);

% Set fracture material properties
fracture.cohesiveLaw     = 'elastic';
fracture.shearStiffness  = 100.0e9;
fracture.normalStiffness = 100.0e9;
fracture.initialAperture = 1.0e-5;
fracture.liquidFluid     = water;

% Add fractures to model
mdl.addPreExistingDiscontinuities(fracture);

%% PROCESS

% Analysis parameters
ti = 0.0;    % Initial time
dt = 1.0;    % Time step
tf = 1.0;   % Final time

% Run analysis
anl = Anl_Transient("Newton");
anl.setUpTransientSolver(ti, dt, tf);
anl.run(mdl);

mdl.plotFieldAlongDiscontinuiy('Aperture',1,'x');
mdl.plotFieldAlongDiscontinuiy('Pressure',1,'x');

% mdl.resetDisplacements();
mdl.resetPressureDirichletBC();
mdl.setUpdateAperture(true);

for i = 1:mdl.nnodes
    if (mdl.NODE(i,1)<1.0e-12 && (abs(mdl.NODE(i,2)-0.5))<0.03)
        mdl.setPressureDirichletBCAtNode(i,11.9e6);
    end
end
mdl.setPressureDirichletBCAtBorder('right', 11.0e6);
mdl.updateDirichletBC();

% clear anl
% anl = Anl_Transient("Newton");
anl.setUpTransientSolver(0.0, 0.01, 500.0, 1, 0.0001, true);

anl.run(mdl);

%% POST-PROCESS

% Plot contours
mdl.plotField('Pressure');

% Plot graphs
Xi = [0.0, 0.0]; Xf = [0.0, 1.0];
% mdl.plotFieldAlongSegment('Pressure', Xi, Xf, 500, 'x');
mdl.plotFieldAlongDiscontinuiy('Aperture',1,'x');
mdl.plotFieldAlongDiscontinuiy('Pressure',1,'x');
mdl.plotFieldAlongDiscontinuiy('Sn',1,'x');
mdl.plotFieldAlongDiscontinuiy('Dn',1,'x');
