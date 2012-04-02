#import the configuration data as read from the GiD
import pfem_var

##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = pfem_var.domain_size


##################################################################
##################################################################
#including kratos path
kratos_benchmarking_path = pfem_var.kratos_path +'/benchmarking' ##kratos_root/benchmarking

import sys
sys.path.append(pfem_var.kratos_path)
sys.path.append(kratos_benchmarking_path)
import benchmarking

##################################################################
##################################################################
 # importing kratos
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.PFEMApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.ExternalSolversApplication import *


##################################################################
##################################################################
def NodeFinder(node_list,X,Y,Z):
   for node in node_list:
	if((node.X-X)**2 + (node.Y-Y)**2 + (node.Z-Z)**2 < .0001):
		return node

def BenchmarkCheck(time, node1):
    benchmarking.Output(time, "Time",0.1,.01)
    print "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
    benchmarking.Output(node1.GetSolutionStepValue(PRESSURE), "Node 1 pressure", 1.0,.01)
    benchmarking.Output(node1.GetSolutionStepValue(VELOCITY_Y), "Node 2 velocity_y", 1.0,.01)


##################################################################
##################################################################
#defining a model part
model_part = ModelPart("FluidPart");  

#importing the solver files and adding the variables
SolverType = pfem_var.SolverType

if(SolverType == "pfem_solver_ale"):
    import pfem_solver_ale
    pfem_solver_ale.AddVariables(model_part)
elif(SolverType == "monolithic_solver_lagrangian"):
    import monolithic_solver_lagrangian
    monolithic_solver_lagrangian.AddVariables(model_part)
else:
    raise "solver type not supported: options are FractionalStep - Monolithic"


#reading a model
name = pfem_var.problem_name

gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO(name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io = ModelPartIO(name)
model_part_io.ReadModelPart(model_part)

##check to ensure that no node has zero density or pressure
for node in model_part.Nodes:
    if(node.GetSolutionStepValue(DENSITY) == 0.0):
        print "node ",node.Id," has zero density!"
        raise 'node with zero density found'
    if(node.GetSolutionStepValue(VISCOSITY) == 0.0):
        print "node ",node.Id," has zero viscosity!"
        raise 'node with zero density found'


mesh_name = 0.0
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh()

print model_part
print model_part.Properties

#setting the limits of the bounding box
box_corner1 = Vector(3);
box_corner1[0]=pfem_var.min_x;
box_corner1[1]=pfem_var.min_y;
box_corner1[2]=pfem_var.min_z;
box_corner2 = Vector(3);
box_corner2[0]=pfem_var.max_x;
box_corner2[1]=pfem_var.max_y;
box_corner2[2]=pfem_var.max_z;

#time setting
output_Dt = pfem_var.output_Dt
max_dt = pfem_var.max_dt
min_dt = pfem_var.min_dt
safety_factor = pfem_var.safety_factor
nsteps = pfem_var.nsteps

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(2)

if(SolverType == "pfem_solver_ale"):
    #adding dofs
    pfem_solver_ale.AddDofs(model_part)
    #creating a fluid solver object
    name = str("dam2d")
    solver = pfem_solver_ale.PFEMSolver(model_part,name,box_corner1,box_corner2,domain_size)
    solver.laplacian_form = pfem_var.laplacian_form
    solver.echo_level = 0
    solver.prediction_order = 1
    solver.predictor_corrector = True
    solver.smooth = True
    solver.alpha_shape = 1.2
    solver.max_press_its = 3;  
    #initializing the solver
    initial_dt = 0.001*min_dt
    solver.Initialize(initial_dt,output_Dt)
elif(SolverType == "monolithic_solver_lagrangian"):
    #adding dofs

    monolithic_solver_lagrangian.AddDofs(model_part)
    solver = monolithic_solver_lagrangian.MonolithicSolver(model_part,domain_size,box_corner1,box_corner2)
    oss_swith = pfem_var.use_oss
    dynamic_tau = pfem_var.dynamic_tau
    solver.echo_level = 2
    model_part.ProcessInfo.SetValue(OSS_SWITCH, oss_swith);				
    model_part.ProcessInfo.SetValue(DYNAMIC_TAU, dynamic_tau);
    solver.Initialize(output_Dt)


###############################################################
back_node = NodeFinder(model_part.Nodes , 0.5137, 0.5526, 0.0)
print back_node

###############################################################
time = 0.0
for step in range(0,nsteps):
    print "line49"

    new_Dt = solver.EstimateDeltaTime(min_dt,max_dt)

    time = time + new_Dt*safety_factor
    
    model_part.CloneTimeStep(time)

    print time

    #solving the fluid problem
    if(step > 3):
        solver.Solve(time,gid_io)
        BenchmarkCheck(time, back_node)


          
        
