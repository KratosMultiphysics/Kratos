%% DESCRIPTION
%
% Uniform traction on a plate with isotropic damage model.
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
mdl.addPorePressure = true;

%% MESH

% Create mesh
Lx = 4000.0;  % Horizontal dimension (m)
Ly = 1000.0;  % Vertical dimension (m)
Nx = 50;      % Number of elements in the x-direction
Ny = 25;      % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny, [], [], 'ISOQ4', 0.5, 0.7, 0.5, 0.0);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Earth pressure coefficient
K0 = 0.65;

% Create porous media
rock = PorousMedia('rock');   
rock.Young = 30.0e+9;        % Young modulus (Pa)
rock.nu    = K0/(1+K0);      % Poisson ratio

% Set materials to model
mdl.setMaterial(rock);

%% BOUNDARY CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('left',   [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('right',  [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('bottom', [NaN, 0.0]);

% Pressure load
mdl.addLoadAtBorder('top', 2, -70.0e6);

% Set external pore-pressure
P = 35.0e6 * ones(mdl.nnodes,1);
mdl.setPorePressureField(P);

%% PROCESS

% Run analysis
anl = Anl_Linear();
anl.run(mdl);

% Update the pressure at the left side
tol = 1.0e-5;
reservoir = isInsideRectangle(mdl.NODE, [0.0-tol,350.0-tol], [0.5*Lx+tol,650.0+tol]);
P(reservoir == 1) = P(reservoir == 1) + 20.0e6;

% Update pressure at the right side
offset = 100.0;
reservoir = isInsideRectangle(mdl.NODE, [0.5*Lx,350.0+offset-tol], [Lx+tol,650.0+offset+tol]);
P(reservoir == 1) = P(reservoir == 1) + 20.0e6;
mdl.setPorePressureField(P);

% Run analysis
anl.run(mdl);

%% POST-PROCESS

% Print stresses
fprintf("\n\n IP      Sxx          Syy          Szz          Sxy\n")
elem = mdl.element(1).type;
for i = 1:elem.nIntPoints
    stress = elem.intPoint(i).stress;
    fprintf("%2d   %+.4e  %+.4e  %+.4e  %+.4e \n",i,stress(1),stress(2),stress(3),stress(4));
end

% Plot stresses
mdl.plotField('Sx');
mdl.plotField('Sy');
mdl.plotField('Sxy');

mdl.plotFieldAlongSegment('Sxy',[0.5*Lx,0.0],[0.5*Lx,Ly],100,'y')
