# Exact solution
def Pressure(node):
  return -2*node.X

def VelocityX(node):
  return node.Y*(2-node.Y)

def VelocityY(node):
  return 0.0

def BenchmarkCheck(VelX,VelY,P1,P2):
    benchmarking.Output(VelX, "x velocity at center node", None, 0.001)
    benchmarking.Output(VelY, "x velocity at center node", None, 0.001)
    benchmarking.Output(P1, "pressure at lower right corner", None, 0.01)
    benchmarking.Output(P2, "pressure at upper left corner", None, 0.01)

##################################################################
kratos_path = '../../../..'
kratos_benchmarking_path    = '../../../../benchmarking' ##kratos_root/benchmarking
import sys
sys.path.append(kratos_path)
sys.path.append(kratos_benchmarking_path)

from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *

import benchmarking

#defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart");  

##################################################################

###
## Problem definition
###

domain_size = 2
problem_name = "trapezoid"
dynamic_tau = 1.0
oss_switch = 1
print_output = False
time_step = 0.1
max_time = 1.0
# assuming Density = 1.0, Viscosity = 1.0
# mdpa assings this, as well as FLAG_VARIABLE=1.0 in the contour

# Import solver and define solution step data
import monolithic_solver_eulerian as monolithic_solver
monolithic_solver.AddVariables(fluid_model_part)

#reading the fluid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.SingleFile
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO(problem_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io_fluid = ModelPartIO(problem_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)

#setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(3)

##adding dofs
monolithic_solver.AddDofs(fluid_model_part)

# Solver initialization    
fluid_solver = monolithic_solver.MonolithicSolver(fluid_model_part,domain_size)
fluid_solver.rel_vel_tol = 1e-5 # default 1e-5
fluid_solver.abs_vel_tol = 1e-7 # default 1e-7
fluid_solver.rel_pres_tol = 1e-5
fluid_solver.abs_pres_tol = 1e-7
fluid_solver.dynamic_tau = dynamic_tau
fluid_solver.oss_switch  = oss_switch
fluid_solver.Initialize()

## Initial and boundary conditions
for node in fluid_model_part.Nodes:
    if node.GetSolutionStepValue(FLAG_VARIABLE) == 1.0:
      node.Fix(VELOCITY_X)
      node.SetSolutionStepValue(VELOCITY_X,0,VelocityX(node))
      node.Fix(VELOCITY_Y)
      if abs(node.X) < 1e-9 and abs(node.Y) < 1e-9:
         node.Fix(PRESSURE)

#settings to be changed
Dt = time_step
final_time = max_time
output_step = 1

if print_output:
  #mesh to be printed
  mesh_name = 0.0
  gid_io.InitializeMesh( mesh_name)
  gid_io.WriteMesh( fluid_model_part.GetMesh() )
  gid_io.FinalizeMesh()

  gid_io.InitializeResults(mesh_name,(fluid_model_part).GetMesh())

time = 0.0
step = 0
out = 0

pressure_node = None
velocity_node = None
for node in fluid_model_part.Nodes:
    if node.Id == 58:
        velocity_node = node
    elif node.Id == 1:
        lower_node = node
    elif node.Id == 100:
        upper_node = node

# initialize solution history (fill solution step data buffer)
while(step < fluid_model_part.GetBufferSize()):
    time = time + Dt
    fluid_model_part.CloneTimeStep(time)
    step = step + 1
t0 = time
final_time += time

while(time < final_time):

    time = time + Dt
    fluid_model_part.CloneTimeStep(time)

    fluid_solver.Solve()

    if print_output and out == output_step:
        gid_io.WriteNodalResults(PRESSURE,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VELOCITY,fluid_model_part.Nodes,time,0)
        out = 0

    out = out + 1
    step = step + 1
 
VelX = velocity_node.GetSolutionStepValue(VELOCITY_X,0)
VelY = velocity_node.GetSolutionStepValue(VELOCITY_Y,0)
P1 = lower_node.GetSolutionStepValue(PRESSURE,0)
P2 = upper_node.GetSolutionStepValue(PRESSURE,0)
BenchmarkCheck(VelX,VelY,P1,P2)

if print_output:
    gid_io.FinalizeResults()

