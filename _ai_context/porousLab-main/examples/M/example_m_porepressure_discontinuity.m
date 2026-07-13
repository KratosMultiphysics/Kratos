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
mdl.condenseEnrDofs   = false;
mdl.addPorePressure   = true;
mdl.subDivIntegration = true;
mdl.symmetricSDAEFEM  = true;

%% MESH

% Create mesh
Lx = 2.0e-1;  % Horizontal dimension (m)
Ly = 2.0e-1;  % Vertical dimension (m)
Nx = 1;       % Number of elements in the x-direction
Ny = 1;       % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create porous media
rock = PorousMedia('rock');   
rock.Young = 30.0e+9;  % Young modulus (Pa)
rock.nu    = 0.2;      % Poisson ratio

% Set materials to model
mdl.setMaterial(rock);

%% BOUNDARY CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('left',   [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('right',  [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('bottom', [NaN, 0.0]);

% Set external pore-pressure
% P = 10 * [0;1;0;1];
P = 10 * [1;0;1;0];
% P = 10 * [1;1;1;1];
mdl.setPorePressureField(P);

%% DISCONTINUITIES

% Create discontinuities 
Xd = [ 0.0 , 0.5*Ly;
       Lx  , 0.5*Ly ];
fault = Discontinuity(Xd, true);

% Set fracture material properties
fault.cohesiveLaw     = 'elastic';
fault.shearStiffness  = 1.0e15;       % Pa/m
fault.normalStiffness = 1.0e15;       % Pa/m

% Add fractures to model
discontinuityData = struct('addTangentialStretchingMode', true, 'addNormalStretchingMode', true, 'addRelRotationMode', false);
mdl.addPreExistingDiscontinuities(fault, discontinuityData);

% Set pressure jump at the discontinuities
nDiscontinuities = mdl.getNumberOfDiscontinuities();
for i = 1:nDiscontinuities
    nDiscontinuitySeg = mdl.discontinuitySet(i).getNumberOfDiscontinuitySegments();
    for j = 1:nDiscontinuitySeg
        mdl.discontinuitySet(i).segment(j).DP = -10;
    end
end

%% PROCESS

% Run analysis
anl = Anl_Linear();
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
mdl.plotFieldAlongDiscontinuiy('Sn',1,'x')
