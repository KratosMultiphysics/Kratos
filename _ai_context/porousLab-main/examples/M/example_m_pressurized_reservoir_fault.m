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
mdl.addPorePressure   = true;

%% MESH

% Create mesh
Lx = 4000.0;  % Horizontal dimension (m)
Ly = 1000.0;  % Vertical dimension (m)
Nx = 60;      % Number of elements in the x-direction
Ny = 60;      % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny, [], [], 'ISOQ4', 0.5, 0.7, 0.5, 0.8);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Earth pressure coefficient
% K0 = 0.60;

% Create porous media
rock = PorousMedia('rock');   
rock.Young = 11868.0e+6;        % Young modulus (Pa)
rock.nu    = 0.29;         % Poisson ratio

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
P0 = 35.0e6;
P = P0 * ones(mdl.nnodes,1);
mdl.setPorePressureField(P);

%% DISCONTINUITIES

% Set formulation options
mdl.symmetricSDAEFEM = false;
mdl.condenseEnrDofs  = false;

% Fault center
Xdc = 0.5 * [Lx, Ly];

% Fault dip
dip = 60.0;
dxdy = 0.0;
if (dip-90)<1.0e-8
    dxdy = 1.0/tand(dip);
end

% Create discontinuities 
Xd = [ Xdc(1)+(0.0 - Xdc(2))*dxdy , 0.0;
       Xdc(1)+(Ly - Xdc(2))*dxdy , Ly ];
fault = Discontinuity(Xd, true);

% Set fracture material properties
fault.cohesiveLaw     = 'elastic';
fault.shearStiffness  = 1.0e15;       % Pa/m
fault.normalStiffness = 1.0e15;       % Pa/m

% Add fractures to model
discontinuityData = struct('addTangentialStretchingMode', false, 'addNormalStretchingMode', false, 'addRelRotationMode', false);
mdl.addPreExistingDiscontinuities(fault);

%% PROCESS

% Run analysis
anl = Anl_Linear();
anl.run(mdl);

mdl.resetDisplacements();

% Parameters
tol    = 1.0e-5;

% Update pressure at the reservoir
reservoir = isInsideRectangle(mdl.NODE, [0.0-tol,350.0-tol], [Lx+tol,650.0+tol]);
P(reservoir == 1) = P0 + 20.0e6;

% Update values
mdl.setPorePressureField(P);

% Run analysis
anl.run(mdl);

%% POST-PROCESS

mdl.plotField('PressureExt');
hold on;
fault.plotIntersectedGeometry();

% Analytical solution of the stresses at the continuum
npts = 50;
K0 = rock.nu/(1.0 - rock.nu);
y_plot = linspace(0.0, Ly, npts);
sy_anal = zeros(npts,1);
for i = 1:npts
    if (y_plot(i) >= 350.0) && (y_plot(i) <= 650.0)
        sy_anal(i) = -70.0e6 + 35.0e6 + 20.0e6; 
    else
        sy_anal(i) = -70.0e6 + 35.0e6;
    end
end
sx_anal = K0 * sy_anal;

% Analytical solution of the stresses along the fault
sn_anal = zeros(npts,1);
st_anal = zeros(npts,1);
s_anal = zeros(npts,1);
dXd = Xd(2,:) - Xd(1,:);
ld = norm(dXd);
md = dXd / ld;
cs = md(1); sn = md(2);
Rot = [cs , sn; -sn , cs];
for i = 1:npts
    coord_d_i = Xd(1,:) + (i - 1) * (ld / (npts - 1)) * md;
    s_anal(i) = (i - 1) * (ld / (npts - 1));
    if (coord_d_i(2) >= 350.0) && (coord_d_i(2) <= 650.0)
        sv = -70.0e6 + 35.0e6 + 20.0e6; 
    else
        sv = -70.0e6 + 35.0e6;
    end
    sh = K0 * sv;
    S = [sh , 0; 0, sv];
    Smn = Rot * S * Rot';
    sn_anal(i) = Smn(2,2);
    st_anal(i) = Smn(1,2);
end

mdl.plotFieldAlongSegment('Sy',[0.5*Lx,0.0],[0.5*Lx,Ly],100,'y'); hold on
plot(sy_anal,y_plot,'ok','DisplayName','Analytical');

mdl.plotFieldAlongSegment('Sx',[0.5*Lx,0.0],[0.5*Lx,Ly],100,'y'); hold on
plot(sx_anal,y_plot,'ok','DisplayName','Analytical');

mdl.plotFieldAlongDiscontinuiy('St',1,'y'); hold on
plot(st_anal,s_anal,'ok','DisplayName','Analytical');

mdl.plotFieldAlongDiscontinuiy('Sn',1,'y'); hold on
plot(sn_anal,s_anal,'ok','DisplayName','Analytical');

mdl.plotField('Uy');
mdl.plotField('Ux');

