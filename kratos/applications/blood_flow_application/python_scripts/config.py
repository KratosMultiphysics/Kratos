#makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division
# PARAMETER_CONFIG
relative_path_3D = "/home/biomedical/Desktop/IamWorkingOn_June/FFR/Test_mesh/"
relative_path_1D = "/home/biomedical/Desktop/IamWorkingOn_June/FFR/1D_Clean/"
name = "test_2"
artery_type=[1]
deactivate_list=[22,23]
inlets_1d=[[100,22]]
outlets_1d=[[1001,23]]
# Config_pressure_conditions
systolic_pressure = 15598.7
diastolic_pressure = 9065.9
#systolic_pressure = 10000
#diastolic_pressure = 10000
time_period = 1
nro_cardiac_cycles = 6
# Config_Hypermia_pressure
systolic_hypermia_pressure = 15598.7
diastolic_hypermia_pressure = 9065.9 #8666.0
#systolic_hypermia_pressure = 10000
#diastolic_hypermia_pressure = 10000 #8666.0
nro_FFR_cardiac_cycles = 3


###makes KratosMultiphysics backward compatible with python 2.6 and 2.7
##from __future__ import print_function, absolute_import, division
### PARAMETER_CONFIG
##relative_path_3D = "/home/biomedical/Desktop/IamWorkingOn_June/FFR/3D/"
##relative_path_1D = "/home/biomedical/Desktop/IamWorkingOn_June/1D_Clean/"
##name = "volume_mesh"
##artery_type = [8]
##deactivate_list=[13,14]
##inlets_1d=[[100,13]]
##outlets_1d=[[1001,14]]
### Config_pressure_conditions
##systolic_pressure = 10000
##diastolic_pressure = 8000
##time_period = 0.02
##nro_cardiac_cycles = 1
### Config_Hypermia_pressure
##systolic_hypermia_pressure = 10000
##diastolic_hypermia_pressure = 5000
##Resistance_factor=0.9
##nro_FFR_cardiac_cycles = 1
