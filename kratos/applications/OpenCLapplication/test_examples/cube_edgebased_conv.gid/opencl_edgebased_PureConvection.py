import time as my_timer
tenter = my_timer.time()

import edgebased_var

##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = edgebased_var.domain_size

from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.OpenCLApplication import *

#including kratos path
kratos_benchmarking_path = '../../../../benchmarking'
import sys
sys.path.append(kratos_benchmarking_path)
import benchmarking

#defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart");  

#############################################


##importing the solvers needed
fluid_model_part.AddNodalSolutionStepVariable(DISTANCE)
##fluid_model_part.AddNodalSolutionStepVariable(TEMP_CONV_PROJ)
fluid_model_part.AddNodalSolutionStepVariable(VELOCITY)
fluid_model_part.AddNodalSolutionStepVariable(AUX_INDEX)

def BenchmarkCheck(time, model_part):
    dist = 0.0; 
    for node in model_part.Nodes:
        if(node.Id == 5420):
            dist = node.GetSolutionStepValue(DISTANCE)
            
    benchmarking.Output(time, "Time")
    benchmarking.Output(dist, "distance on node #5420", 0.00001)



#introducing input file name
input_file_name = edgebased_var.problem_name

#reading the fluid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.SingleFile
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteConditions

##selecting output format
gid_io = GidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
    
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)
print "model part io"

#setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(2)

#neighbour search
number_of_avg_elems = 10
number_of_avg_nodes = 10
neighbour_search = FindNodalNeighboursProcess(fluid_model_part,number_of_avg_elems,number_of_avg_nodes)
(neighbour_search).Execute()

#############################################
##perform refinement
refinement_steps = 2
Refine = LocalRefineTetrahedraMesh(fluid_model_part)
for i in range(0,refinement_steps):
    for elem in fluid_model_part.Elements:
        elem.SetValue(SPLIT_ELEMENT,True)
    refine_on_reference = False;
    interpolate_internal_variables = False;
    Refine.LocalRefineMesh(refine_on_reference,interpolate_internal_variables)    
    (neighbour_search).Execute()

##adding dofs
for node in fluid_model_part.Nodes:
    node.AddDof(DISTANCE);
    node.AddDof(VELOCITY_X);
    node.AddDof(VELOCITY_Y);
    node.AddDof(VELOCITY_Z);
    
print "adddof"
vel = Vector(3);
for node in fluid_model_part.Nodes:
    vel[0] = -node.Y
    vel[1] = node.X
    vel[2] = 0.0
    node.SetSolutionStepValue(VELOCITY,0,vel);

print "assigned vel"
#assigning a cone shaped distance distribution
xc = 1.00/8.00
yc = 1.00/8.00
zc = 0.0
import math
for node in fluid_model_part.Nodes:
    X1 = (node.X - xc) 
    X2 = (node.Y - yc) 
    X3 = (node.Z - zc)
    dist = math.sqrt((X1**2 + X2**2 +  X3**2)) - 0.1
##    if(dist > 0.0):
##        dist = 0.0
    node.SetSolutionStepValue(DISTANCE,0,dist)

t1 = my_timer.time()
    
print "assigned dist dirts"
#constructing the solver
single_device_flag = True
device_group = OpenCLDeviceGroup(cl_device_type.CL_DEVICE_TYPE_GPU,single_device_flag)
device_group.AddCLSearchPath('../../custom_utilities')

matrix_container = OpenCLMatrixContainer3D(device_group)
matrix_container.ConstructCSRVector(fluid_model_part)
matrix_container.BuildCSRData(fluid_model_part)
print "contructed matrix_containers"

##matrix_container.Clear()
##matrix_container = OpenCLMatrixContainer3D(device_group)
##matrix_container.ConstructCSRVector(fluid_model_part)
##matrix_container.BuildCSRData(fluid_model_part)

convection_solver = OpenCLPureConvectionEdgeBased3D(matrix_container, fluid_model_part)
print "solver"
convection_solver.Initialize()
print "initizlized"

t2 = my_timer.time()

#settings to be changed
max_Dt = edgebased_var.max_time_step 
initial_Dt = 0.001 * max_Dt 
final_time = edgebased_var.max_time
output_dt = edgebased_var.output_dt
safety_factor = edgebased_var.safety_factor

number_of_inital_steps = edgebased_var.number_of_inital_steps
initial_time_step = edgebased_var.initial_time_step
out = 0

ProcessInfo = fluid_model_part.ProcessInfo
ProcessInfo
###mesh to be printed
mesh_name = 0.0
gid_io.InitializeMesh( mesh_name)
gid_io.WriteMesh( fluid_model_part.GetMesh() )
gid_io.FinalizeMesh()
gid_io.Flush()

gid_io.InitializeResults(mesh_name, (fluid_model_part).GetMesh());

max_safety_factor = safety_factor
    
time = 0.0
step = 0
next_output_time = output_dt
while(time < final_time):

    convection_solver.ComputeTimeStep(safety_factor)
    Dt = ProcessInfo[DELTA_TIME]
    #print Dt
    time = time + Dt
    fluid_model_part.CloneTimeStep(time)

    #print "******** CURRENT TIME = ",time

    if(step >= 3):
        ##convect levelset function
        convection_solver.Solve()
        
        ## Uncomment to get benchmarking data
##        BenchmarkCheck(time, fluid_model_part)
        
        if(time >=  next_output_time):
            
            gid_io.WriteNodalResults(VELOCITY,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(DISTANCE,fluid_model_part.Nodes,time,0)
            gid_io.Flush()

            next_output_time = time + output_dt

            out = 0

    out = out + 1
    step = step + 1
gid_io.Flush()      
gid_io.FinalizeResults()    
        
t3 = my_timer.time()
print "initial time (reading, refinining etc) =", t1-tenter
print "solver setup time                      =", t2-t1
print "run time                               =", t3-t2
