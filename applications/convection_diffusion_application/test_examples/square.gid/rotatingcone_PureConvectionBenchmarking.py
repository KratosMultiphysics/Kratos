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
    max_temp = 0.0; 
    id_max_temp = 0

    for node in model_part.Nodes:
        temp = node.GetSolutionStepValue(TEMPERATURE);
        if(temp > max_temp):
            max_temp = temp
            id_max_temp = node.Id;

            
        
    benchmarking.Output(time, "Time")
    benchmarking.Output(max_temp, "minimum temperature", 0.00001)
    benchmarking.Output(id_max_temp, "Id of the node with maximum temperature", 0.0)


#defining a model part
model_part = ModelPart("FluidPart");

##importing the solver files and adding the variables
import pure_convection_solver
pure_convection_solver.AddVariables(model_part)
model_part.AddNodalSolutionStepVariable(TEMPERATURE)
model_part.AddNodalSolutionStepVariable(TEMP_CONV_PROJ)
model_part.AddNodalSolutionStepVariable(VELOCITY)
model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
model_part.AddNodalSolutionStepVariable(NODAL_AREA)


##...aqui lista variables para utilizar

#adding of Variables to Model Part should be here when the "very fix container will be ready"
 
#reading a model
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO("square",gid_mode,multifile,deformed_mesh_flag, write_conditions)
gid_io.ReadModelPart(model_part)

print model_part

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(3)

##add Degrees of Freedom to all of the nodes
for node in model_part.Nodes:
    node.AddDof(TEMPERATURE)
    
#settings to be changed
#INITIALIZING FLUID
#assigning the fluid properties

#assigning a rotational velocity field
vel = Vector(3);
for node in model_part.Nodes:
    vel[0] = -node.Y
    vel[1] = node.X
    node.SetSolutionStepValue(VELOCITY,0,vel);

#assigning a cone shaped distance distribution
xc = 1.00/6.00
yc = 1.00/6.00
sigma = 0.2
import math
for node in model_part.Nodes:
    X1 = (node.X - xc) / sigma
    X2 = (node.Y - yc) / sigma
    if( (X1**2 + X2**2) <= 1.00 ):
        temp= 0.25*(1.00+math.cos(math.pi*X1) )*(1.00+math.cos(math.pi*X2) )
        node.SetSolutionStepValue(TEMPERATURE,0,temp)

 


#settings to be changed
output_step = 20
delta_t = 2.00*math.pi/200.0;
out = 0


time_old_print = 0
time = 0.0
max_time = 7.0
step = 0


##########    convection solver
convection_order = 2 #order of the time scheme of the convection solver
pConvPrecond = DiagonalPreconditioner()
convection_linear_solver = linear_solver =  BICGSTABSolver(1e-9, 5000,pConvPrecond)
##convection_solver = PureConvectionUtilities2D();
convection_solver = pure_convection_solver.PureConvectionSolver(model_part,domain_size)

##computing the neighbours
neighbour_finder = FindNodalNeighboursProcess(model_part,10,10);
neighbour_finder.Execute(); ##at wish ... when it is needed

##Variable to be convected:
## #1 = TEPERATURE
## #2 = DISTANCE
convection_solver.scalar_var_convected = 1

convection_solver.Initialize();
print "pure convection solver initialized"

while time < max_time:
    
    time = step*delta_t
    print "Current time = ",time
    model_part.CloneTimeStep(time)

   
    if(step>3):     
        convection_solver.Solve(TEMPERATURE)
        if (benchmarking.InBuildReferenceMode()):
            BenchmarkCheck(time, model_part)
        else:
            BenchmarkCheck(time, model_part)

            
    #print the results
    mesh_name = time #if we want the mesh to change at each time step then ****mesh_name = time****

    step = step + 1

    time_to_print = time - time_old_print

    if(time_to_print >= 0.5 ):
        print "output"
        gid_io.InitializeMesh( mesh_name)
        gid_io.WriteMesh( model_part.GetMesh() )
        gid_io.FinalizeMesh()
        gid_io.InitializeResults(mesh_name , model_part.GetMesh())
##        gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(TEMPERATURE,model_part.Nodes,time,0)

        gid_io.Flush()
        gid_io.FinalizeResults()

        
        time_old_print = time



print "finito"

