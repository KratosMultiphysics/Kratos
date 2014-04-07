from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#
#
# setting the domain size for the problem to be solved
domain_size = 2

#
#
# ATTENTION: here the order is important

# including kratos path
kratos_benchmarking_path = '../../../../benchmarking'  # kratos_root/benchmarking
import sys
sys.path.append(kratos_benchmarking_path)


from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.GPUSolversApplication import *
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


# defining a model part
print("before creation of the model part")
model_part = ModelPart("FluidPart");
print("after creation of the model part")

# importing the solver files and adding the variables
import incompressible_fluid_solver
incompressible_fluid_solver.AddVariables(model_part)

# adding of Variables to Model Part should be here when the "very fix container will be ready"

# reading a model
gid_mode = GiDPostMode.GiD_PostBinary
use_multifile = MultiFileFlag.MultipleFiles
deformed_print_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteConditions
gid_io = GidIO("cavity2d", gid_mode, use_multifile, deformed_print_flag, write_conditions)
write_conditions = WriteConditionsFlag.WriteElementsOnly
# gid_io.ReadMesh(model_part.GetMesh())
gid_io.ReadModelPart(model_part)
# gid_io.WriteMesh((model_part).GetMesh(),domain_size,0.0,GiDPostMode.GiD_PostBinary);
gid_io.InitializeMesh(0.0);
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh()
print(model_part)

# the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(3)

# add Degrees of Freedom to all of the nodes
incompressible_fluid_solver.AddDofs(model_part)


# creating a fluid solver object
fluid_solver = incompressible_fluid_solver.IncompressibleFluidSolver(model_part, domain_size)
fluid_solver.laplacian_form = 3;
fluid_solver.predictor_corrector = True
fluid_solver.vel_toll = 1e-3
fluid_solver.time_order = 2
fluid_solver.echo_level = 2


# pILUPrecond = ILU0Preconditioner()
# fluid_solver.pressure_linear_solver =  BICGSTABSolver(1e-9, 5000,pILUPrecond)
# pDiagPrecond = DiagonalPreconditioner()
# fluid_solver.velocity_linear_solver =  BICGSTABSolver(1e-9, 5000,pDiagPrecond)
# fluid_solver.pressure_linear_solver =  BICGSTABSolver(1e-9, 5000,pDiagPrecond)
# fluid_solver.velocity_linear_solver = SkylineLUFactorizationSolver();
# fluid_solver.pressure_linear_solver = SkylineLUFactorizationSolver();
preSweeps = Vector(10)
postSweeps = Vector(10)
preSweeps[0] = 1
preSweeps[1] = 2
preSweeps[2] = 2
preSweeps[3] = 2
preSweeps[4] = 2
preSweeps[5] = 2
preSweeps[6] = 2
preSweeps[7] = 2
preSweeps[8] = 2
preSweeps[9] = 2

postSweeps[0] = 1
postSweeps[1] = 2
postSweeps[2] = 2
postSweeps[3] = 2
postSweeps[4] = 2
postSweeps[5] = 2
postSweeps[6] = 2
postSweeps[7] = 2
postSweeps[8] = 2
postSweeps[9] = 2
is_preconditioner = True

precond = KratosAMGPreconditioner(4.0 / 3.0, 3, True, 10, 1000, preSweeps, postSweeps, is_preconditioner)
# precond = AMGPreconditioner()
fluid_solver.velocity_linear_solver = GPUBICGSTABSolver(1e-9, 5000)
fluid_solver.pressure_linear_solver = GPUCGSolver(1e-9, 5000, precond)
# fluid_solver.velocity_linear_solver =  GPUBICGSTABSolverWithDiagonalPreconditioner(1e-9, 5000)
# fluid_solver.pressure_linear_solver =  GPUBICGSTABSolverWithDiagonalPreconditioner(1e-9, 5000)
fluid_solver.Initialize()

# settings to be changed
Re = 100.0
nsteps = 100
output_step = 1

Dt = 1000.0 / Re
if(Dt > 0.1):
    Dt = 0.1
out = 0
for node in model_part.Nodes:
    node.SetSolutionStepValue(VISCOSITY, 0, 1.0 / Re)

for node in model_part.Nodes:
    if(node.X < 0.001 and node.Y < 0.001):
        node.Fix(PRESSURE)


zero = Vector(3);
zero[0] = 0.0;
zero[1] = 0.0;
zero[2] = 0.0;

gid_io.InitializeResults(0.0, model_part.GetMesh())

for step in range(1, nsteps):

    time = Dt * step
    print(time)
    model_part.CloneTimeStep(time)

    print(time)
    # print model_part.ProcessInfo()[TIME]

    # solving the fluid problem
    if(step > 3):
        fluid_solver.Solve()
        if (benchmarking.InBuildReferenceMode()):
            BenchmarkCheck(time, model_part)
        else:
            BenchmarkCheck(time, model_part)

    # print the results
    if(out == output_step):
        gid_io.WriteNodalResults(PRESSURE, model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(VELOCITY, model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(PRESS_PROJ, model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(CONV_PROJ, model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(VISCOSITY, model_part.Nodes, time, 0)
        out = 0
    out = out + 1

node = model_part.Nodes[1]

gid_io.FinalizeResults();

print(node)
