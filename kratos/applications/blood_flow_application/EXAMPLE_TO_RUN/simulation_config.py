# SIMULATION_CONFIG
#
# Config_Blood_Properties
blood_density = 1060.0
blood_static_viscosity = 0.0035
blood_viscosity = 0.0035 / 1060.0
#
# Config_Save Results
ascii_results = False
save_results = 0.0001   # save_results = 0.0 (save all results)
#
# Config_Initial_conditions
inlet_pressure_type = "coseno"  # coseno, parabolic,table
Q_initial = 0.0000001
#
# Config_initial_radius
FitRadius = False
#
# Catheter
Use_Catheter = True
Catheter_Radius = 0.00015
#
# Config_Time_step
Sub_steping = True
sub_step = 50  # config.sub_step
step_size_control = True  # True-->Activate	False-->fix (step_size)
step_size = 0.00001  # config.step_size
CardiacCycleConvergence = False
#
# Coupling conditions
FitValues = True  # To get A-B of the 3D model
Coupled_Simulation = True
#
# Config_factor_variables_to_check_different_options(TBR)
pressure_factor = 1
flow_factor = 1
Resistence_factor = 1
Condition_Variable = True
only1Dtest = False  # Use A-B parameters setted in Simulation_config
nodeAB = 0  # Force 1 node A-B Parameters
A = 245630339.384
B = 7.83737812517e+14
