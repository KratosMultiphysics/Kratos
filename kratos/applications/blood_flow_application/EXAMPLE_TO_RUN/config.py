#makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division
# PARAMETER_CONFIG
#relative_path_3D = "/home/biomedical/Desktop/IamWorkingOn_April/Example_TORUN/3D/"
relative_path_3D = "/home/esoudah/kratos/applications/blood_flow_application/EXAMPLE_TO_RUN/3D/"
#relative_path_1D = "/home/biomedical/Desktop/IamWorkingOn_April/Example_TORUN/1D/"
relative_path_1D = "/home/esoudah/kratos/applications/blood_flow_application/EXAMPLE_TO_RUN/1D/"
name = "volume_mesh_tube"
artery_type = [3]
deactivate_list = [40, 41]
inlets_1d = [[100, 40]]
outlets_1d = [[1001, 41]]
# Config_pressure_conditions
systolic_pressure = 12500.0
diastolic_pressure = 7500.0
time_period = 0.1
nro_cardiac_cycles = 2
# Config_Hypermia_pressure
diastolic_hypermia_pressure = 7500.0
