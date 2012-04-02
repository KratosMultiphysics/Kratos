def BenchmarkCheck(Vel303x,Vel303y,Vel312,Vel322):
    benchmarking.Output(Vel303x, "x velocity at node 303", None, 1e-5)
    benchmarking.Output(Vel303y, "y velocity at node 303", None, 1e-5)
    benchmarking.Output(Vel312, "x velocity at node 312", None, 1e-5)
    benchmarking.Output(Vel322, "x velocity at node 312", None, 1e-5)

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
problem_name = "slip_test"
dynamic_tau = 1.0
oss_switch = 0
print_output = False
time_step = 0.1
max_time = 1.0
y_wall = 0.01
density = 1.0
viscosity = 0.01

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
fluid_model_part.SetBufferSize(2)

##adding dofs
monolithic_solver.AddDofs(fluid_model_part)

# Solver initialization    
fluid_solver = monolithic_solver.MonolithicSolver(fluid_model_part,domain_size)
rel_vel_tol = 1e-5 # default 1e-5
abs_vel_tol = 1e-7 # default 1e-7
rel_pres_tol = 1e-5
abs_pres_tol = 1e-7

# The solver is initialized here (instead of calling fluid_solver.Initialize()
# to use a non-default time scheme
fluid_solver.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent( fluid_solver.alpha,fluid_solver.move_mesh_strategy,domain_size )
fluid_solver.conv_criteria = VelPrCriteria(rel_vel_tol,abs_vel_tol,rel_pres_tol,abs_pres_tol)
fluid_solver.solver = ResidualBasedNewtonRaphsonStrategy(fluid_solver.model_part,\
                                                         fluid_solver.time_scheme,\
                                                         fluid_solver.linear_solver,
                                                         fluid_solver.conv_criteria,
                                                         fluid_solver.max_iter,
                                                         fluid_solver.CalculateReactionFlag,
                                                         fluid_solver.ReformDofSetAtEachStep,
                                                         fluid_solver.MoveMeshFlag)
(fluid_solver.solver).SetEchoLevel(fluid_solver.echo_level)

fluid_model_part.ProcessInfo.SetValue(DYNAMIC_TAU, dynamic_tau);
fluid_model_part.ProcessInfo.SetValue(OSS_SWITCH, oss_switch );
fluid_model_part.ProcessInfo.SetValue(M, fluid_solver.regularization_coef );

# Calculate normals
normal_calculator = NormalCalculationUtils()
normal_calculator.CalculateOnSimplex(fluid_model_part,domain_size)

# Assing IS_STRUCTURE to slip nodes
for condition in fluid_model_part.Conditions:
  if condition.GetValue(IS_STRUCTURE) == 1.0:
    for node in condition.GetNodes():
      node.SetValue(IS_STRUCTURE,1.0)
      node.SetValue(Y_WALL,y_wall)
      #node.SetSolutionStepValue(IS_STRUCTURE,0,1.0) # This is not used by Kratos, added here to view it in the results.
    condition.SetValue(IS_STRUCTURE,0.0)

#result nodes
node303 = None
node312 = None
node322 = None

for node in fluid_model_part.Nodes:
  node.SetSolutionStepValue(DENSITY,0,density)
  node.SetSolutionStepValue(VISCOSITY,0,viscosity)
  # "improve" normals on the corners
  if node.Id == 1:# or node.Id == 480:
    node.SetSolutionStepValue(NORMAL_X,0,1.0)
    node.SetSolutionStepValue(NORMAL_Y,0,0.0)
    node.SetSolutionStepValue(NORMAL_Z,0,0.0)
  elif node.Id == 103:# or node.Id == 704:
    node.SetSolutionStepValue(NORMAL_X,0,-1.0)
    node.SetSolutionStepValue(NORMAL_Y,0,0.0)
    node.SetSolutionStepValue(NORMAL_Z,0,0.0)
  elif node.Id == 480 or node.Id == 704:
    node.SetSolutionStepValue(VELOCITY_X,0,0.0)
    node.SetSolutionStepValue(VELOCITY_Y,0,0.0)
    node.SetSolutionStepValue(VELOCITY_Z,0,0.0)
  elif node.Id == 303:
    node303 = node
  elif node.Id == 312:
    node312 = node
  elif node.Id == 322:
    node322 = node

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
        gid_io.WriteNodalResults(NORMAL,fluid_model_part.Nodes,time,0)
#        gid_io.WriteNodalResults(VELOCITY,fluid_model_part.Nodes,time,0)
        out = 0

    out = out + 1
    step = step + 1
 
Vel303x = node303.GetSolutionStepValue(VELOCITY_X,0)
Vel303y = node303.GetSolutionStepValue(VELOCITY_Y,0)
Vel312 = node312.GetSolutionStepValue(VELOCITY_X,0)
Vel322 = node322.GetSolutionStepValue(VELOCITY_X,0)
BenchmarkCheck(Vel303x,Vel303y,Vel312,Vel322)

if print_output:
    gid_io.FinalizeResults()

