##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
##### copyright by CIMNE, Barcelona, Spain                   #####
#####          and Janosch Stascheit for TUNCONSTRUCT        #####
##### all rights reserved                                    #####
##################################################################
#setting the domain size for the problem to be solved
domain_size = 3
##################################################################
##################################################################
## ATTENTION: here the order is important                    #####
##################################################################
## including kratos path                                     #####
## ATTENTION: the following lines have to be adapted to      #####
##            match your acrtual configuration               #####
##################################################################
import sys
import os
kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##setting up paths
kratos_libs_path = kratos_root_path+'libs' ##kratos_root/libs
kratos_applications_path = kratos_root_path+'applications' ##kratos_root/applications
##################################################################
##################################################################
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

##################################################################
##################################################################
sys.path.append('./consolidation_test.gid')
import consolidation_test_include
from consolidation_test_include import *
# calculate insitu-stress for geology_virgin.gid
model = consolidation_test_include.Model('consolidation_test','./')
linear_elastic = False
model.InitializeModel( linear_elastic )

##################################################################
###  SIMULATION  #################################################
##################################################################
# =====================
# | USER SCRIPT FOR CALCULATION OF AUSBLAS.GID |
# vvvvvvvvvvvvvvvvvvvvv
#simulation script for consolidation test

#loading phase
pressure = 0.0
delta_pressure = -1000.0
for node in model.layer_nodes_sets['surface']:
    model.model_part.Nodes[node].SetSolutionStepValue( FACE_LOAD_Z, pressure )
time = 0.0
delta_time = 0.0001
for step in range(0,10):
    pressure = pressure + delta_pressure
    for node in model.layer_nodes_sets['surface']:
        model.model_part.Nodes[node].SetSolutionStepValue( FACE_LOAD_Z, pressure )
    time = time + delta_time
    model.Solve( time, 0, 0, 0, 0 )
    model.WriteOutput( time )
#consolidation phase

for node in model.layer_nodes_sets['surface']:
    #model.model_part.Nodes[node].SetSolutionStepValue( FACE_LOAD_Z, pressure )
    model.model_part.Nodes[node].Fix( WATER_PRESSURE )
    model.model_part.Nodes[node].SetSolutionStepValue( WATER_PRESSURE_NULL, 0.0 )
    model.model_part.Nodes[node].SetSolutionStepValue( WATER_PRESSURE_EINS, 0.0 )
    model.model_part.Nodes[node].SetSolutionStepValue( WATER_PRESSURE, 0.0 )
for step in range(0,10):
    time = time + delta_time
    model.Solve( time, 0, 0, 0, 0 )
    model.WriteOutput( time )
delta_time = 0.001
for step in range(0,10):
    time = time + delta_time
    model.Solve( time, 0, 0, 0, 0 )
    model.WriteOutput( time )
delta_time = 0.01
for step in range(0,10):
    time = time + delta_time
    model.Solve( time, 0, 0, 0, 0 )
    model.WriteOutput( time )
delta_time = 0.1
for step in range(0,10):
    time = time + delta_time
    model.Solve( time, 0, 0, 0, 0 )
    model.WriteOutput( time )
delta_time = 0.5
for step in range(0,50):
    time = time + delta_time
    model.Solve( time, 0, 0, 0, 0 )
    model.WriteOutput( time )
delta_time = 1.0
for step in range(0,50):
    time = time + delta_time
    model.Solve( time, 0, 0, 0, 0 )
    model.WriteOutput( time )
delta_time = 5.0
for step in range(0,50):
    time = time + delta_time
    model.Solve( time, 0, 0, 0, 0 )
    model.WriteOutput( time )
print("all done.")
