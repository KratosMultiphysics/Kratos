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

import benchmarking

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

#importing applications
import applications_interface
applications_interface.Import_ConvectionDiffusionApplication = True
applications_interface.Import_IncompressibleFluidApplication = True
applications_interface.Import_ExternalSolversApplication = False
applications_interface.ImportApplications(kernel, kratos_applications_path)

## from now on the order is not anymore crucial
##################################################################
##################################################################
from KratosConvectionDiffusionApplication import *
from KratosIncompressibleFluidApplication import *
##from KratosR1ExternalSolversApplication import *

##################################################################
##################################################################
def NodeFinder(node_list,X,Y,Z):
   for node in node_list:
	if((node.X-X)**2 + (node.Y-Y)**2 + (node.Z-Z)**2 < .000001):
		return node


def BenchmarkCheck(time, node1):
    benchmarking.Output(time, "Time")
    benchmarking.Output(node1.GetSolutionStepValue(PRESSURE), "Node 1 pressure", 0.05)
    benchmarking.Output(node1.GetSolutionStepValue(DISTANCE), "Node 1 distance", 0.05)


#defining a model part
model_part = ModelPart("FluidPart");

##importing the solver files and adding the variables
import level_set_elembased_fluid_solver
level_set_elembased_fluid_solver.AddVariables(model_part)
##...aqui lista variables para utilizar

#adding of Variables to Model Part should be here when the "very fix container will be ready"

#reading a model
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = EdgebasedGidIO("QuietWater_LevelSet",gid_mode,multifile,deformed_mesh_flag, write_conditions)
gid_io.ReadModelPart(model_part)

print model_part

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(3)

##add Degrees of Freedom to all of the nodes
level_set_elembased_fluid_solver.AddDofs(model_part)

#settings to be changed
#INITIALIZING FLUID
zero = Vector(3);
zero[0] = 0.0;  zero[1] = 0.0;  zero[2] = 0.0;
body_force = Vector(3);
body_force[0] = 0.0;    body_force[1] = -9.806; body_force[2] = 0.0;
ext_press = 0.0

velocity = zero[0] 
pressure = ext_press

###############################################################
back_node = NodeFinder(model_part.Nodes , 4.0, 2.0, 0.0)
print back_node
###############################################################

for node in model_part.Nodes:

    node.SetSolutionStepValue(DENSITY,0,1000)
    node.SetSolutionStepValue(VISCOSITY,0,1e-6)
##    node.SetSolutionStepValue(VELOCITY_X,0,0.0)
##    node.SetSolutionStepValue(VELOCITY_Y,0,0.0)
##    node.SetSolutionStepValue(VELOCITY_Z,0,0.0)
    node.SetSolutionStepValue(BODY_FORCE_Y,0,body_force[1])
##      node.SetSolutionStepValue(SEEPAGE_DRAG, 0,zero)
    node.Free(PRESSURE)

    node.SetSolutionStepValue(PRESSURE,0,pressure)
    node.SetSolutionStepValue(EXTERNAL_PRESSURE,0,0.0)

##necessary setting of the distances values as IC
for node in model_part.Nodes:
    node.SetSolutionStepValue(DISTANCE,0,node.Y-5.0)
    node.SetSolutionStepValue(DIAMETER,0,1.0)
    node.SetSolutionStepValue(POROSITY,0,1.0)
    
delta_t = 0.05

time_old_print = 0.0
time = 0.0
max_time = 10.0
step = 0


extrapolation_distance = 1
#FluidSolver****************************************************************
fluid_solver = level_set_elembased_fluid_solver.ElemBasedLevelSetSolver(model_part,domain_size,body_force)
#PressureSolver*************************************************************
pPrecond = DiagonalPreconditioner()
linear_solver = BICGSTABSolver(1e-3,5000,pPrecond)
##linear_solver = CGSolver(1e-5,5000,pPrecond)
##linear_solver = SkylineLUFactorizationSolver()
##linear_solver = SuperLUSolver()

fluid_solver.redistance_frequency  = 1
fluid_solver.predictor_corrector = True
fluid_solver.reform_convection_matrix = True
fluid_solver.ReformDofAtEachIteration = True



fluid_solver.number_of_extrapolation_layers = 3

#computing the neighbours
neighbour_finder = FindNodalNeighboursProcess(model_part,10,10);
neighbour_finder.Execute(); ##at wish ... when it is needed

fluid_solver.Initialize()
print "fluid solver initialized"

fluid_solver.model_part.ProcessInfo.SetValue(DYNAMIC_TAU,1)

while time < max_time:
    
    time = step*delta_t    
    print "Current time = ",time
    model_part.CloneTimeStep(time)

    if(step > 3):
        
        fluid_solver.Solve()
        BenchmarkCheck(time, back_node)
        #print the results
    mesh_name = time #if we want the mesh to change at each time step then ****mesh_name = time****

    step = step + 1

    time_to_print = time - time_old_print

    if(time_to_print >= 0.0 ):
        print "output"
        gid_io.InitializeMesh( mesh_name)
        gid_io.WriteMesh( model_part.GetMesh() )
        gid_io.FinalizeMesh()
        gid_io.InitializeResults(mesh_name , model_part.GetMesh())
        gid_io.WriteNodalResults(PRESSURE,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DISTANCE,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DENSITY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DIAMETER, model_part.Nodes,time,0)
##        gid_io.WriteNodalResults(NORMAL,model_part.Nodes,time,0)
##        gid_io.WriteNodalResults(POROSITY,model_part.Nodes,time,0)
####        gid_io.WriteNodalResults(IS_BOUNDARY,model_part.Nodes,time,0)
##        gid_io.WriteNodalResults(IS_STRUCTURE,model_part.Nodes,time,0)
####        gid_io.WriteNodalResults(IS_FLUID,model_part.Nodes,time,0)
####        gid_io.WriteNodalResults(SEEPAGE_DRAG, model_part.Nodes, time, 0)
##        gid_io.WriteNodalResults(CONVECTION_VELOCITY, model_part.Nodes, time, 0)
####        gid_io.WriteNodalResults(TEMPERATURE,model_part.Nodes,time,0)
##        gid_io.WriteNodalResults(VISCOSITY,model_part.Nodes,time,0)
##        gid_io.WriteNodalResults(BODY_FORCE,model_part.Nodes,time,0)
##        gid_io.WriteNodalResults(IS_DIVIDED,model_part.Nodes,time,0)


        gid_io.Flush()
        gid_io.FinalizeResults()

        
        time_old_print = time



print "finito"

