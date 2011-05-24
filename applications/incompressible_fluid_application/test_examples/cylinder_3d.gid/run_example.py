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
kratos_benchmarking_path = fluid_only_var.kratos_path +'/benchmarking' ##kratos_root/benchmarking
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
applications_interface.Import_IncompressibleFluidApplication = True
applications_interface.Import_ExternalSolversApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)

## from now on the order is not anymore crucial
##################################################################
##################################################################
from KratosIncompressibleFluidApplication import *
from KratosExternalSolversApplication import *

##################################################################
##################################################################
def NodeFinder(node_list,X,Y,Z):
   for node in node_list:
	if((node.X-X)**2 + (node.Y-Y)**2 + (node.Z-Z)**2 < .001):
		return node

def BenchmarkCheck(time, node1):
    benchmarking.Output(time, "Time")
    benchmarking.Output(node1.GetSolutionStepValue(PRESSURE), "Node 1 pressure", 1.0)
    benchmarking.Output(node1.GetSolutionStepValue(VELOCITY_Y), "Node 2 velocity_y", 1.0)


#defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart");  

#############################################
##importing the solvers needed
SolverType = fluid_only_var.SolverType
if(SolverType == "FractionalStep"):
    import fractional_step_solver
    fractional_step_solver.AddVariables(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian"):
    import monolithic_solver_eulerian
    monolithic_solver_eulerian.AddVariables(fluid_model_part)
else:
    raise "solver type not supported: options are FractionalStep - Monolithic"

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
fluid_model_part.SetBufferSize(3)

##adding dofs
if(SolverType == "FractionalStep"):
    fractional_step_solver.AddDofs(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian"):
    monolithic_solver_eulerian.AddDofs(fluid_model_part)


#########select here the laplacian form!!!!!!!!!!!!!!!!!
laplacian_form = fluid_only_var.laplacian_form 
if(laplacian_form >= 2):
    for node in fluid_model_part.Nodes:
        node.Free(PRESSURE)

##check to ensure that no node has zero density or pressure
for node in fluid_model_part.Nodes:
    if(node.GetSolutionStepValue(DENSITY) == 0.0):
        print "node ",node.Id," has zero density!"
        raise 'node with zero density found'
    if(node.GetSolutionStepValue(VISCOSITY) == 0.0):
        print "node ",node.Id," has zero viscosity!"
        raise 'node with zero VISCOSITY found'

#creating the solvers
#fluid solver
if(SolverType == "FractionalStep"):
    fluid_solver = fractional_step_solver.IncompressibleFluidSolver(fluid_model_part,domain_size)
    fluid_solver.laplacian_form = laplacian_form; #standard laplacian form
    fluid_solver.predictor_corrector = fluid_only_var.predictor_corrector
    fluid_solver.max_press_its = fluid_only_var.max_press_its
    fluid_solver.Initialize()
elif(SolverType == "monolithic_solver_eulerian"): 
    fluid_solver = monolithic_solver_eulerian.MonolithicSolver(fluid_model_part,domain_size)
    oss_swith = fluid_only_var.use_oss
    dynamic_tau = fluid_only_var.dynamic_tau
    fluid_model_part.ProcessInfo.SetValue(OSS_SWITCH, oss_swith);				
    fluid_model_part.ProcessInfo.SetValue(DYNAMIC_TAU, dynamic_tau);
    fluid_solver.Initialize()


print "fluid solver created"

#settings to be changed
Dt = fluid_only_var.Dt 
full_Dt = Dt 
initial_Dt = 1.0* full_Dt #0.05 #0.01
final_time = fluid_only_var.max_time
output_step = fluid_only_var.output_step

out = 0


#mesh to be printed
mesh_name = 0.0
gid_io.InitializeMesh( mesh_name)
gid_io.WriteMesh( fluid_model_part.GetMesh() )
gid_io.FinalizeMesh()

gid_io.InitializeResults(mesh_name,(fluid_model_part).GetMesh())

###############################################################
back_node = NodeFinder(fluid_model_part.Nodes , 8.04104, 2.30717, 1.18733)
print back_node
print "maaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaammmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm"

###############################################################


time = 0.0
step = 0
while(time < final_time):

    if(step < 5):
        Dt = initial_Dt
    else:
        Dt = full_Dt
        
    time = time + Dt
    fluid_model_part.CloneTimeStep(time)

    if(step >= 3):
        fluid_solver.Solve()
        BenchmarkCheck(time, back_node)

    if(out == output_step):
        gid_io.WriteNodalResults(PRESSURE,fluid_model_part.Nodes,time,0)
	gid_io.WriteNodalResults(AIR_PRESSURE,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VELOCITY,fluid_model_part.Nodes,time,0)
        #gid_io.WriteNodalResults(DISPLACEMENT,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(MESH_VELOCITY,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(IS_STRUCTURE,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(IS_BOUNDARY,fluid_model_part.Nodes,time,0)
	gid_io.WriteNodalResults(IS_POROUS,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(IS_FREE_SURFACE,fluid_model_part.Nodes,time,0)
        #gid_io.PrintOnGaussPoints(THAWONE,fluid_model_part,time)
        #gid_io.PrintOnGaussPoints(THAWTWO,fluid_model_part,time)
        gid_io.WriteNodalResults(ADVPROJ,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DIVPROJ,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DENSITY,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DENSITY_AIR,fluid_model_part.Nodes,time,0)
       ## gid_io.WriteNodalResults(NODAL_H,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VISCOSITY,fluid_model_part.Nodes,time,0)

        out = 0

    out = out + 1
    step = step + 1
      
gid_io.FinalizeResults()
          
        

