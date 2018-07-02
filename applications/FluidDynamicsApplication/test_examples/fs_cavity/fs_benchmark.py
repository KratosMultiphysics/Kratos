from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
problem_name = "cavity2D"
Dt = 0.1
final_time = 10.0
write_post_file = False  # True
output_step = 1

#
#
# setting the domain size for the problem to be solved
domain_size = 2

import sys
sys.path.append('../../../../')
sys.path.append('../../../../benchmarking')
from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *

import benchmarking


def BenchmarkCheck(time, model_part):
    max_press = 0.0
    min_press = 0.0
    vel2min = 10000.0
    id_min_vel = 0
    x_min_vel = 0.0
    y_min_vel = 0.0
    for node in model_part.Nodes:
        press = node.GetSolutionStepValue(PRESSURE)
        if(press > max_press):
            max_press = press
        elif(press < min_press):
            min_press = press

        x = node.X
        y = node.Y
        vel = node.GetSolutionStepValue(VELOCITY)
        vel2 = vel[0] ** 2 + vel[1] ** 2
        if(x > 0.1 and x < 0.9 and y > 0.1 and y < 0.9):
            if(vel2 < vel2min):
                vel2min = vel2
                id_min_vel = node.Id
                x_min_vel = node.X
                y_min_vel = node.Y

    benchmarking.Output(time, "Time",1e-7)
    benchmarking.Output(min_press, "minimum pressure", 0.00001)
    benchmarking.Output(max_press, "maximum pressure", 0.00001)
    benchmarking.Output(id_min_vel, "Id of the node with minimum velocity norm", 0.0)
    benchmarking.Output(x_min_vel, "coord x minimum velocity norm", 0.0)
    benchmarking.Output(y_min_vel, "coord y minimum velocity norm", 0.0)

# defining a model part for the fluid and one for the structure
model = Model()
fluid_model_part = model.CreateModelPart("FluidPart")

#
# importing the solvers needed
import vms_fractional_step_solver as fractional_step_solver
fractional_step_solver.AddVariables(fluid_model_part)

# reading the fluid part
model_part_io_fluid = ModelPartIO(problem_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(3)

# adding dofs
fractional_step_solver.AddDofs(fluid_model_part)
# check to ensure that no node has zero density or pressure
for node in fluid_model_part.Nodes:
    if(node.GetSolutionStepValue(DENSITY) == 0.0):
        print("node ", node.Id, " has zero density!")
        raise 'node with zero density found'
    if(node.GetSolutionStepValue(VISCOSITY) == 0.0):
        print("node ", node.Id, " has zero viscosity!")
        raise 'node with zero VISCOSITY found'

# creating the solvers
# fluid solver
fluid_solver = fractional_step_solver.IncompressibleFluidSolver(fluid_model_part, domain_size)
fluid_solver.predictor_corrector = False
fluid_solver.ReformDofAtEachIteration = False
fluid_solver.max_press_its = 3
fluid_solver.max_vel_its = 5
fluid_solver.velocity_linear_solver = SkylineLUFactorizationSolver();
fluid_solver.pressure_linear_solver = SkylineLUFactorizationSolver();
fluid_solver.Initialize()
print("fluid solver created")

# mesh to be printed
if write_post_file:
    gid_mode = GiDPostMode.GiD_PostBinary
    multifile = MultiFileFlag.MultipleFiles
    deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
    write_conditions = WriteConditionsFlag.WriteElementsOnly
    gid_io = GidIO(problem_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)

    mesh_name = 0.0
    gid_io.InitializeMesh(mesh_name)
    gid_io.WriteMesh(fluid_model_part.GetMesh())
    gid_io.FinalizeMesh()

    gid_io.InitializeResults(mesh_name, (fluid_model_part).GetMesh())


out = 0
time = 0.0
step = 0
while(time < final_time):

    time = time + Dt
    fluid_model_part.CloneTimeStep(time)

    if(step >= 3):
        fluid_solver.Solve()
        BenchmarkCheck(time, fluid_model_part)

    if write_post_file and out == output_step:
        gid_io.WriteNodalResults(PRESSURE, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(VELOCITY, fluid_model_part.Nodes, time, 0)
        out = 0

    out = out + 1
    step = step + 1

if write_post_file:
    gid_io.FinalizeResults()
