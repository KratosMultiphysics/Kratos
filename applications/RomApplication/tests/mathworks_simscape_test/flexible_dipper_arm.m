function arm = flexible_dipper_arm()
clear all
close all
clc
%% Step 1: Define the Geometry and Material Properties of the Dipper Arm
% Geometry
stlFile = 'sm_flexible_dipper_arm_bin.stl';
stl_obj = stlread(stlFile);
% figure(1)
% trisurf(stl_obj);
% axis equal

% Material
E = 200e9;     % Young's modulus in Pa
nu = 0.26;     % Poisson's ratio (nondimensional)
rho = 7800;    % Mass density in kg/m^3

%% Step 2: Specify the Locations of Interface Frames
origins = [-0.500  0      0      % Frame 1: Cylinder connection point
            1.500  0      0      % Frame 2: Bucket connection point
            0     -0.130  0];    % Frame 3: Fulcrum point
numFrames = size(origins,1);

%% Step 3: Create the Finite-Element Mesh
feModel = createpde('structural','modal-solid');% returns a structural analysis model of , 'modal-solid' to create a structural model for modal analysis of a solid (3-D) problem.
importGeometry(feModel,stlFile); %importGeometry(model,geometryfile) creates a geometry container from the specified STL geometry file, and includes the geometry in the model container.
structuralProperties(feModel, ...
    'YoungsModulus',E, ...
    'PoissonsRatio',nu, ...
    'MassDensity',rho);
generateMesh(feModel, ...
    'GeometricOrder','quadratic', ...
    'Hmax',0.2, ...
    'Hmin',0.02);
size(feModel.Mesh.Nodes);
%% Step 4: Set up the Multipoint Constraints for the Interface Frames
% Identify the geometric regions 
% figure
% pdegplot(feModel,'FaceLabels','on','FaceAlpha',0.5)

faceIDs = [1,27,23];    % List in the same order as the interface frame origins

% To verify these values, plot the mesh and highlight the selected faces:
% figure
% pdemesh(feModel,'FaceAlpha',0.5)
% hold on
% colors = ['rgb' repmat('k',1,numFrames-3)];
% assert(numel(faceIDs) == numFrames);
% for k = 1:numFrames
%     nodeIdxs = findNodes(feModel.Mesh,'region','Face',faceIDs(k));
%     scatter3( ...
%         feModel.Mesh.Nodes(1,nodeIdxs), ...
%         feModel.Mesh.Nodes(2,nodeIdxs), ...
%         feModel.Mesh.Nodes(3,nodeIdxs), ...
%         'ok','MarkerFaceColor',colors(k))
%     scatter3( ...
%         origins(k,1), ...
%         origins(k,2), ...
%         origins(k,3), ...
%         80,colors(k),'filled','s')
% end
% hold off

% Call the function structuralBC (Partial Differential Equation Toolbox) to define the MPCs for the boundary nodes in these faces:
for k = 1:numFrames
    structuralBC(feModel,'Face',faceIDs(k),'Constraint','multipoint', 'Reference', origins(k,:));
end

%% Step 5: Generate the Reduced-Order Model
% Craig-Bampton order reduction method
rom = reduce(feModel,'FrequencyRange',[0,1e3]);
% Store the results of the reduction in a data structure arm. 
arm.P = rom.ReferenceLocations';  % Interface frame locations (n x 3 matrix)
arm.K = rom.K;                    % Reduced stiffness matrix
arm.M = rom.M;                    % Reduced mass matrix

% Compute modal damping Matrix
dampingRatio = 0.05;
arm.C = computeModalDampingMatrix(dampingRatio,rom.K,rom.M);

% Reorder as the corresponding interface frames on the block
frmPerm = zeros(numFrames,1);    % Frame permutation vector
dofPerm = 1:size(arm.K,1);       % DOF permutation vector

assert(size(arm.P,1) == numFrames);
for i = 1:numFrames
    for j = 1:numFrames
        if isequal(arm.P(j,:),origins(i,:))
            frmPerm(i) = j;
            dofPerm(6*(i-1)+(1:6)) = 6*(j-1)+(1:6);
            continue;
        end
    end
end
assert(numel(frmPerm) == numFrames);
assert(numel(dofPerm) == size(arm.K,1));

arm.P = arm.P(frmPerm,:);
arm.K = arm.K(dofPerm,:);
arm.K = arm.K(:,dofPerm);
arm.M = arm.M(dofPerm,:);
arm.M = arm.M(:,dofPerm);
arm.C = arm.C(dofPerm,:);
arm.C = arm.C(:,dofPerm);
%% Save stiffness matrix for python to compare with Kratos'
% filename="matriz%d.mat";
% str = compose(filename,Nodos(2));
% str = 'matriz.mat';
% matriz=arm.K;
% save(str,'matriz')
end