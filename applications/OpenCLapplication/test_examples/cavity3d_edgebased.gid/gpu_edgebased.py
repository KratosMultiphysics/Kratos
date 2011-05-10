# -*- coding: utf-8 -*-
import time as my_timer
tenter = my_timer.time()

import fluid_only_var

##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = fluid_only_var.domain_size

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path            = fluid_only_var.kratos_path + '/libs' ##kratos_root/libs
kratos_applications_path    = fluid_only_var.kratos_path + '/applications' ##kratos_root/applications
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

#importing applications
import applications_interface
applications_interface.Import_IncompressibleFluidApplication = True
applications_interface.Import_OpenCLApplication = True
applications_interface.Import_MeshingApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)

## from now on the order is not anymore crucial
##################################################################
##################################################################
from KratosIncompressibleFluidApplication import *
from KratosOpenCLApplication import *
from KratosMeshingApplication import *
##from KratosExternalSolversApplication import *

#defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart");  

#############################################
##importing the solvers needed
import opencl_eulerian_NS_solver
opencl_eulerian_NS_solver.AddVariables(fluid_model_part)
fluid_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)

#introducing input file name
input_file_name = fluid_only_var.problem_name

#reading the fluid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)


#setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(2)

#neighbour search
number_of_avg_elems = 10
number_of_avg_nodes = 10
neighbour_search = FindNodalNeighboursProcess(fluid_model_part,number_of_avg_elems,number_of_avg_nodes)
(neighbour_search).Execute()

#############################################
##perform refinement
refinement_steps = 1
Refine = LocalRefineTetrahedraMesh(fluid_model_part)
for i in range(0,refinement_steps):
    for elem in fluid_model_part.Elements:
        elem.SetValue(SPLIT_ELEMENT,True)
    refine_on_reference = False;
    interpolate_internal_variables = False;
    Refine.LocalRefineMesh(refine_on_reference,interpolate_internal_variables)    
    (neighbour_search).Execute()

##aux_renumberer = RenumberByNeighbourCountUtil()
##aux_renumberer.Renumber(fluid_model_part.Nodes)

opencl_eulerian_NS_solver.AddDofs(fluid_model_part)

t1 = my_timer.time()

zero = Vector(3)
zero[0]=0.0;
zero[1]=0.0;
zero[2]=0.0;
density = 1.0
viscosity = 0.01
safety_factor = 0.9

#creating the solvers
#fluid solver
cl_sources='../../custom_utilities'
fluid_solver = opencl_eulerian_NS_solver.OpenClSolver(fluid_model_part,domain_size,zero,viscosity,density,cl_sources)
fluid_solver.Initialize()

print "fluid solver created"

#settings to be changed
Dt = fluid_only_var.Dt 
full_Dt = Dt 
initial_Dt = 0.001 * full_Dt #0.05 #0.01
final_time = fluid_only_var.max_time
output_step = fluid_only_var.output_step

out = 0

#mesh to be printed
mesh_name = 0.0
gid_io.InitializeMesh( mesh_name)
gid_io.WriteMesh( fluid_model_part.GetMesh() )
gid_io.FinalizeMesh()

gid_io.InitializeResults(mesh_name,(fluid_model_part).GetMesh())

t2 = my_timer.time()

time = 0.0
step = 0
while(time < final_time):

    if(step < 3):
        Dt = initial_Dt
    else:
        Dt = fluid_solver.EstimateTimeStep(safety_factor,full_Dt)

    print "current time = ",time, " Dt = ",Dt
    time = time + Dt
    fluid_model_part.CloneTimeStep(time)

    if(step >= 3):
        fluid_solver.Solve()

    if(out == output_step):
        gid_io.WriteNodalResults(PRESSURE,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VELOCITY,fluid_model_part.Nodes,time,0)

        out = 0

    out = out + 1
    step = step + 1
      
gid_io.FinalizeResults()
    
t3 = my_timer.time()
print "initial time (reading, refinining etc) =", t1-tenter
print "solver setup time                      =", t2-t1
print "run time                               =", t3-t2

        

