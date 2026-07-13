%% DESCRIPTION
%
% Strip footing problem.
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
mdl.intOrder = 2;  % Integration quadrature order

%% MESH

% Load mesh
load('MeshStripFooting');
[node, elem] = convertToQuadraticMesh(node, elem);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create porous media
rock = PorousMedia('rock');
rock.mechanical    = 'druckerPrager';  % Constitutive law
rock.MCmatch       = 'planestrain';    % How Drucker-Prager surfaces matches Mohr-Coulomb
rock.Young         = 1.0e+7;           % Young modulus (kPa)
rock.nu            = 0.48;             % Poisson ratio
rock.cohesion      = 490.0;            % Cohesion (kPa)
rock.frictionAngle = 20*pi/180;        % Friction angle (rad)
rock.dilationAngle = 20*pi/180;        % Dilation angle (rad)

% Set materials to model
mdl.setMaterial(rock);

%% BOUNDARY CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('left',   [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('right',  [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('bottom', [NaN, 0.0]);

for i = 1:mdl.nnodes
    if ((mdl.NODE(i,1) <= 0.2) && (mdl.NODE(i,2) == 5.0))
        mdl.addLoadAtNode(i, [0.0 , -1.0e2])
    end
end

%% PROCESS

% Setup analysis
anl = Anl_NonlinearQuasiStatic('GeneralizedDisplacement');
anl.adjustStep    = true;
anl.increment     = 0.01;
anl.max_increment = 0.1;
anl.max_lratio    = 10.0;
anl.max_step      = 10;
anl.max_iter      = 100;
anl.trg_iter      = 3;

% Node and DOF used to plot Load Factor vs Displacement
ndId = mdl.closestNodeToPoint([0.5, 5.0]);
anl.setPlotDof(ndId, 2);

% Run analysis
anl.echo = false;
anl.run(mdl);