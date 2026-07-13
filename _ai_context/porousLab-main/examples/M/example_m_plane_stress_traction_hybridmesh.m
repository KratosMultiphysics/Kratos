%% DESCRIPTION
%
% Uniform traction on a plate with linear elastic model.
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
node = [0.0 , 0.0;
        0.0 , 1.0;
        1.0 , 0.0;
        1.0 , 1.0;
        1.5 , 0.5;
        2.0 , 0.0;
        2.0 , 1.0];

elem = {[1, 3, 4, 2];
        [3, 6, 5];
        [6, 7, 5];
        [7, 4, 5];
        [4, 3, 5]};

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create porous media
rock = PorousMedia('rock');
rock.mechanical = 'elastic';  % Constitutive law
rock.Young      = 2.0e+10;    % Young modulus (Pa)
rock.nu         = 0.0;        % Poisson ratio

% Set materials to model
mdl.setMaterial(rock);

%% BOUNDARY CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('left', [0.0, 0.0]);

% Loads
mdl.addLoadAtBorder('right', 1, 2.0e+6);

%% PROCESS

% Run analysis
anl = Anl_Linear();
anl.run(mdl);

%% POS-PROCESS

% Print results to command window
mdl.printResults();

% Plot contours
mdl.plotFieldAlongSegment('Ux',[0.0,0.5],[2.0,0.5],10)
mdl.plotField('Ux');
mdl.plotField('Sx');
