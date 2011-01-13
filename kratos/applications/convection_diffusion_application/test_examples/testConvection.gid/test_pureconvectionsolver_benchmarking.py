##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 2

import math
##import cProfile
##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path = '../../../../libs' ##kratos_root/libs
kratos_applications_path = '../../../../applications/' ##kratos_root/applications
kratos_benchmarking_path = '../../../../benchmarking' ##kratos_root/benchmarking
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)
sys.path.append(kratos_benchmarking_path)

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

#importing applications
import applications_interface
applications_interface.Import_ConvectionDiffusionApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)

## from now on the order is not anymore crucial
##################################################################
##################################################################
from KratosConvectionDiffusionApplication import *
import benchmarking

def BenchmarkCheck(time, model_part):
    dist = 0.0; 
    for node in model_part.Nodes:
        if(node.Id == 15):
            dist = node.GetSolutionStepValue(DISTANCE)
            
    benchmarking.Output(time, "Time")
    benchmarking.Output(dist, "distance on node #15", 0.00001)


#defining a model part
model_part = ModelPart("FluidPart");

##########################################################
thermal_settings = ConvectionDiffusionSettings()
thermal_settings.SetDensityVariable(DENSITY)
thermal_settings.SetDiffusionVariable(CONDUCTIVITY)
thermal_settings.SetUnknownVariable(TEMPERATURE)
thermal_settings.SetVolumeSourceVariable(HEAT_FLUX)
thermal_settings.SetSurfaceSourceVariable(FACE_HEAT_FLUX)
thermal_settings.SetMeshVelocityVariable(MESH_VELOCITY)
##########################################################

import pure_convection_solver
pure_convection_solver.AddVariables(model_part,thermal_settings)

##importing the solver files and adding the variables
import pure_convection_solver
pure_convection_solver.AddVariables(model_part,thermal_settings)
model_part.AddNodalSolutionStepVariable(VELOCITY)
model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
model_part.AddNodalSolutionStepVariable(TEMP_CONV_PROJ)
model_part.AddNodalSolutionStepVariable(DISTANCE)
model_part.AddNodalSolutionStepVariable(NODAL_AREA)
##...aqui lista variables para utilizar

#adding of Variables to Model Part should be here when the "very fix container will be ready"
 
#reading a model
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO("testConvection",gid_mode,multifile,deformed_mesh_flag, write_conditions)
gid_io.ReadModelPart(model_part)

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(3)

#importing the solver files
pure_convection_solver.AddDofs(model_part,thermal_settings)

##add DEGREES OF FREEDOM to all of the nodes: here you can choose the variable to be convected!!
for node in model_part.Nodes:
    node.AddDof(DISTANCE)


    
#settings to be changed
#INITIALIZING FLUID
zero = Vector(3);
zero[0] = 0.0;  zero[1] = 0.0;  zero[2] = 0.0;
ext_press = 0.0

velocity = zero[0] 
pressure = ext_press

for node in model_part.Nodes:
    node.SetSolutionStepValue(VELOCITY_X,0,1.0)
    node.SetSolutionStepValue(VELOCITY_Y,0,0.0)
    node.SetSolutionStepValue(VELOCITY_Z,0,0.0)
    node.SetSolutionStepValue(DISTANCE,0,node.X)

delta_t = 0.002

time_old_print = 0
time = 0.0
max_time = 5.0
step = 0



##########    convection solver: CONSTRUCTOR ##########
convection_order = 2
pConvPrecond = DiagonalPreconditioner()
convection_linear_solver = linear_solver =  BICGSTABSolver(1e-9, 5000,pConvPrecond)
##convection_solver = PureConvectionUtilities2D();
convection_solver = pure_convection_solver.PureConvectionSolver(model_part,domain_size,thermal_settings)

##computing the neighbours
neighbour_finder = FindNodalNeighboursProcess(model_part,10,10);
neighbour_finder.Execute(); ##at wish ... when it is needed
##

convection_solver.number_of_extrapolation_layers = 3
convection_solver.ReformDofAtEachIteration = True
##Variable to be convected:
## #1 = TEPERATURE
## #2 = DISTANCE
convection_solver.scalar_var_convected = 2



##########    convection solver: INITIALIZE ##########
convection_solver.Initialize();



while time < max_time:
    
    time = step*delta_t
    print "Current time = ",time
    model_part.CloneTimeStep(time)


    if(step>=3):
        
        convection_solver.Solve(DISTANCE)
        if (benchmarking.InBuildReferenceMode()):
            BenchmarkCheck(time, model_part)
        else:
            BenchmarkCheck(time, model_part)
      
    #print the results
    mesh_name = time #if we want the mesh to change at each time step then ****mesh_name = time****

    step = step + 1

    time_to_print = time - time_old_print

    if(time_to_print >= 0.2 ):
        print "output"
        gid_io.InitializeMesh( mesh_name)
        gid_io.WriteMesh( model_part.GetMesh() )
        gid_io.FinalizeMesh()
        gid_io.InitializeResults(mesh_name , model_part.GetMesh())
        gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DISTANCE,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(TEMP_CONV_PROJ, model_part.Nodes,time,0)
        gid_io.WriteNodalResults(MESH_VELOCITY,model_part.Nodes,time,0)



        gid_io.Flush()
        gid_io.FinalizeResults()

        
        time_old_print = time



print "finito"

