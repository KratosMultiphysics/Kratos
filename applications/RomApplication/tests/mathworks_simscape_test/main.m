close all; clear all; clc;
tic
%% Path to python-Kratos
% pathToKratos = fileparts(which("MainKratos.py"));
% if count(py.sys.path,pathToKratos) == 0
%     insert(py.sys.path,int32(0),pathToKratos);
% end
% system("python3 MainKratos.py");
%% Name of simulink project
model_name = "flexible_dipper_arm_kratos";
%% Load block diagram flexible dipper arm
flexible_dipper_arm_kratos % Name of simscape/simulink project
%% Run Flexible dipper arm analysis
arm = flexible_dipper_arm(); % To obtain mass matrix and damping matrix, not yet implemented in Kratos
%% Run Kratos simulation
fname = 'Stiffness_Matrix.json'; % Read stiffness matrix json file created by Kratos' simulation
K = jsondecode(fileread(fname)); % Save stiffness matrix on current workspace
arm.K = K.StiffnessMatrix; % Replace Simscape's stiffness matrix with Kratos' 
arm.P = K.InterfaceFrameOrigins; % Replace Simscape's interface frame origins with Kratos'
%% Modify Simscape model
mdlWks = get_param(model_name,'ModelWorkspace'); % Save Simscape's Model Workspace
mdlWks.setVariablePart('arm.K',arm.K); % Replace Stiffness Matrix
mdlWks.setVariablePart('arm.P',arm.P); % Replace Interface Frame Origins
%% Set the translation of the Rigid Transform blocks appropriately
P12 = arm.P(1,:) - arm.P(2,:); % Compute the relative offset between the interface frames and the common
P13 = arm.P(1,:) - arm.P(3,:); % reference frame

% Use the following command to set the value of the Cartesian offset
% translation
set_param(model_name+'/Rigid Transform F1-F2','TranslationCartesianOffset', ['[' num2str(P12) ']']);
set_param(model_name+'/Rigid Transform F1-F3','TranslationCartesianOffset', ['[' num2str(P13) ']']);
%% Run Simscape model with new parameters 
sim(model_name);
%% Plot the Simscape's nodes results
figure
title('Displacements Frame 1')
xlabel('Time (s)')
ylabel('Displacement (mm)')
subplot(3,1,1)
plot(DXF2); hold on; plot(DYF2); hold on; plot(DZF2)
legend('X Disp','Y Disp','Z Disp')
title('Displacements Frame 2')
xlabel('Time (s)')
ylabel('Displacement (mm)')
subplot(3,1,2)
plot(DXF3); hold on; plot(DYF3); hold on; plot(DZF3)
legend('X Disp','Y Disp','Z Disp')
title('Displacements Frame 3')
xlabel('Time (s)')
ylabel('Displacement (mm)')
subplot(3,1,3)
plot(QF2); hold on; plot(QF3)
legend('Frame 2','Frame 3')
title('Frame Quaternions')
xlabel('Time (s)')
ylabel('Displacement (mm)')
%% Export data to Kratos
sim_res.Time = DXF2.time; % Obtain time steps array
sim_res.Num_of_frames = length(arm.P);
sim_res.Kratos_Interface_Frame_Origins = K.InterfaceFrameOrigins;
sim_res.Kratos_Displacement = K.Displacement; % Displacement of Kratos' process
sim_res.Kratos_Rotation = K.Rotation; % Rotation of Kratos' process in radians
for i=2:sim_res.Num_of_frames
    temp_name = strcat( 'Interface_',num2str(i)); % Dynamic Field-Names
    sim_res.(temp_name).Disp_x = eval(strcat('DXF',num2str(i))).data; % Obtain displacements on x array
    sim_res.(temp_name).Disp_y = eval(strcat('DYF',num2str(i))).data; % Obtain displacements on y array
    sim_res.(temp_name).Disp_z = eval(strcat('DZF',num2str(i))).data; % Obtain displacements on z array
    sim_res.(temp_name).Quaternion = eval(strcat('QF',num2str(i))).data; % Obtain quaternions array
end
save('sim_res.mat','-struct','sim_res') % Create a .mat file with the results structure
toc
%% Run Kratos' post-proces HROM
% system('python3 MainKratos_HROM.py');