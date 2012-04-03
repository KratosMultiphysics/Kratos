##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 2

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_path = '../../../../' ##kratos_root/
kratos_benchmarking_path = '../../../../benchmarking' ##kratos_root/benchmarking

import sys
sys.path.append(kratos_lpath)
sys.path.append(kratos_benchmarking_path)

#importing Kratos main library
from KratosMultiphysics import *

#importing applications
from KratosMultiphysics.IncompressibleFluidApplicationimport *
from KratosMultiphysics.PFEMApplication import *
from KratosMultiphysics.MeshingApplication import *

import benchmarking

## from now on the order is not anymore crucial
##################################################################
##################################################################

def FindNode(node_list,x,y,z):
    for node in node_list:
        if ((node.X - x) ** 2 + (node.Y - y) ** 2 + (node.Z - z) ** 2 < 0.0000001):
            return node
    
def BenchmarkCheck(time, node1, node2):
    benchmarking.Output(time, "Time")
    benchmarking.Output(node1.GetSolutionStepValue(PRESSURE), "Node 1 pressure", 1.0)
    benchmarking.Output(node2.GetSolutionStepValue(PRESSURE), "Node 2 pressure", 1.0)


#defining a model part
model_part = ModelPart("FluidPart");
print "aaa"

#importing the solver files and adding the variables
import pfem_solver_ale
pfem_solver_ale.AddVariables(model_part)

#reading a model
gid_mode_flag = GiDPostMode.GiD_PostBinary
use_multifile = MultiFileFlag.MultipleFiles
deformed_print_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteConditions
gid_io = GidIO("dam2d",gid_mode_flag,use_multifile,deformed_print_flag,write_conditions)
gid_io.ReadModelPart(model_part)
print model_part

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(3)

#adding dofs
pfem_solver_ale.AddDofs(model_part)

#setting the limits of the bounding box
box_corner1 = Vector(3); box_corner1[0]=-0.1; box_corner1[1]=-0.1; box_corner1[2]=-0.1;
box_corner2 = Vector(3); box_corner2[0]=1.1; box_corner2[1]=1.1;  box_corner2[2]=-0.0;

#creating a fluid solver object
name = str("dam2d")
solver = pfem_solver_ale.PFEMSolver(model_part,name,box_corner1,box_corner2,domain_size)
solver.laplacian_form = 2
solver.echo_level = 0
solver.prediction_order = 1
solver.predictor_corrector = True
solver.smooth = True
solver.alpha_shape = 1.2
solver.max_press_its = 3;  
#solver.vel_toll = 1e-12
#solver.max_vel_its = 20;
##solver.echo_level = 2
solver.pressure_linear_solver = SkylineLUFactorizationSolver()
solver.velocity_linear_solver = SkylineLUFactorizationSolver()


gravity = Vector(3)
gravity[0] = 0.00; gravity[1] = -9.81; gravity[2] = 0.0;
for node in model_part.Nodes:
    node.SetSolutionStepValue(VISCOSITY,0,0.000001)
    node.SetSolutionStepValue(DENSITY,0,1000.00000)
    node.SetSolutionStepValue(BODY_FORCE,0,gravity)
 
Dt = 0.01
nsteps = 1000
output_Dt = 0.05
min_dt = 0.005
max_dt = 0.02

safety_factor = 0.5;

##########################

node_1 = FindNode(model_part.Nodes, 0.5, 0.0, 0.0)
node_2 = FindNode(model_part.Nodes, 0.24, 0.0, 0.0)

##########################

#initializing the solver
solver.Initialize(Dt,output_Dt)

time = Dt
for step in range(0,nsteps):
    print "solution step =" , step

    new_Dt = solver.EstimateDeltaTime(min_dt,max_dt)
   # print "time = ", time, " new_Dt= ",new_Dt," step = ", step

    #new_Dt = Dt

    time = time + new_Dt*safety_factor

    #time = Dt*step
    model_part.CloneTimeStep(time)

    print time
    #print model_part.ProcessInfo()[TIME]

    #solving the fluid problem
    print "qui"
    if(step > 3):
        solver.Solve(time,gid_io)
##        if(step > 4):
##            solver.box_corner2[1] = 0.1
    print "li"
    
    BenchmarkCheck(time, node_1, node_2)

print "solution finished"



          
        

