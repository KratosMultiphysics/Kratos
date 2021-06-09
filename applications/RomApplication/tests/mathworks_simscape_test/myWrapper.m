function time_step = myWrapper(varargin)
%% Modify Simscape model
fname = 'Stiffness_Matrix.json'; % Read stiffness matrix json file created by Kratos' simulation
K = jsondecode(fileread(fname)); % Save stiffness matrix on current workspace
K_k = K.StiffnessMatrix;%*(1+varargin{1}); % Replace Simscape's stiffness matrix with Kratos'
mdlWks = get_param("my_sm_flexible_dipper_arm",'ModelWorkspace'); % Save Simscape's Model Workspace
mdlWks.setVariablePart('arm.K',K_k); % Replace Stiffness Matrix
Num_of_frames = size(K.InterfaceFrameOrigins,1);
%% Export data to Kratos
%mod = py.importlib.reload(py.importlib.import_module("MainKratos_HROM"));
%mod.RunKratosHROM(varargin{:},Num_of_frames);
py.MainKratos_HROM.RunKratosHROM(varargin{:},Num_of_frames);
%% Copy new postprocess into folder 
new = sprintf('Kratos_postprocess/VISUALIZE_HROM_0_%d.vtk',varargin{1});
movefile('vtk_output/VISUALIZE_HROM_0_1.vtk',new);
time_step = varargin{1}+1;
disp(time_step)