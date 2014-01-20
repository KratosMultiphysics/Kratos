from __future__ import unicode_literals, print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# import Kratos
import ProjectParameters
import sys
import os
from ctypes import *
sys.path.append(ProjectParameters.kratos_path)
from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.EmpireApplication import *

# EMPIRE_API & connect
print('subSystem 2')
libempire_api = cdll.LoadLibrary(os.environ['EMPIRE_API_LIBSO_ON_MACHINE'])
libempire_api.EMPIRE_API_Connect("subSystem2.xml")

domain_size = ProjectParameters.domain_size

# defining a model part for the fluid
fluid_model_part = ModelPart("FluidPart")

# add variables to model part
fluid_model_part.AddNodalSolutionStepVariable(REACTION)
fluid_model_part.AddNodalSolutionStepVariable(DISTANCE)
fluid_model_part.AddNodalSolutionStepVariable(PRESSURE)
fluid_model_part.AddNodalSolutionStepVariable(IS_INTERFACE)

# importing the solver needed
import vms_fractional_step_solver as solver
solver.AddVariables(fluid_model_part)

# introducing input file name
input_file_name = ProjectParameters.problem_name

# reading the fluid part

# initialize GiD  I/O
gid_mode = GiDPostMode.GiD_PostBinary

if ProjectParameters.GiDWriteMeshFlag:
    deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
else:
    deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed

if(ProjectParameters.VolumeOutput):
    if ProjectParameters.GiDWriteConditionsFlag:
        write_conditions = WriteConditionsFlag.WriteConditions
    else:
        write_conditions = WriteConditionsFlag.WriteElementsOnly
else:
    write_conditions = WriteConditionsFlag.WriteConditions

multifile = MultiFileFlag.SingleFile

gid_io = GidIO(input_file_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(3)

# adding dofs
solver.AddDofs(fluid_model_part)

# If Lalplacian form = 2, free all pressure Dofs
laplacian_form = ProjectParameters.laplacian_form
if(laplacian_form >= 2):
    for node in fluid_model_part.Nodes:
        node.Free(PRESSURE)

# decoupled solver schemes: need solver for pressure and solver fol velocity

if ProjectParameters.SolverType in ["FractionalStep"]:
    # Velocity preconditioner
    try:
        if ProjectParameters.Velocity_Preconditioner_type == 'Diagonal':
            vel_precond = DiagonalPreconditioner()
        elif ProjectParameters.Velocity_Preconditioner_type == 'ILU0':
            vel_precond = ILU0Preconditioner()
    except AttributeError:
        # If we are using a direct solver Velocity_Preconditioner_type will not
        # exist. This is not an error.
        vel_precond = None

    # Velocity solver
    if ProjectParameters.Velocity_Linear_Solver == "Conjugate gradient":
        if ProjectParameters.Velocity_Preconditioner_type == 'None':
            velocity_linear_solver = CGSolver(ProjectParameters.Velocity_Iterative_Tolerance, ProjectParameters.Velocity_Solver_Max_Iteration)
        else:
            velocity_linear_solver = CGSolver(ProjectParameters.Velocity_Iterative_Tolerance, ProjectParameters.Velocity_Solver_Max_Iteration, vel_precond)
    elif ProjectParameters.Velocity_Linear_Solver == "BiConjugate gradient stabilized":
        if ProjectParameters.Velocity_Preconditioner_type == 'None':
            velocity_linear_solver = BICGSTABSolver(ProjectParameters.Velocity_Iterative_Tolerance, ProjectParameters.Velocity_Solver_Max_Iteration)
        else:
            velocity_linear_solver = BICGSTABSolver(ProjectParameters.Velocity_Iterative_Tolerance, ProjectParameters.Velocity_Solver_Max_Iteration, vel_precond)
    elif ProjectParameters.Velocity_Linear_Solver == "GMRES":
        if ProjectParameters.Velocity_Preconditioner_type == 'None':
            velocity_linear_solver = GMRESSolver(ProjectParameters.Velocity_Iterative_Tolerance, ProjectParameters.Velocity_Solver_Max_Iteration)
        else:
            velocity_linear_solver = GMRESSolver(ProjectParameters.Velocity_Iterative_Tolerance, ProjectParameters.Velocity_Solver_Max_Iteration, vel_precond)
    elif ProjectParameters.Velocity_Linear_Solver == "Skyline LU factorization":
        velocity_linear_solver = SkylineLUFactorizationSolver()
    elif ProjectParameters.Velocity_Linear_Solver == "Super LU":
        velocity_linear_solver = SuperLUSolver()
    elif ProjectParameters.Velocity_Linear_Solver == "Parallel MKL Pardiso":
        from KratosMultiphysics.MKLSolversApplication import MKLPardisoSolver
        velocity_linear_solver = MKLPardisoSolver()

    # Pressure preconditioner
    try:
        if ProjectParameters.Pressure_Preconditioner_type == 'Diagonal':
            press_precond = DiagonalPreconditioner()
        elif ProjectParameters.Pressure_Preconditioner_type == 'ILU0':
            press_precond = ILU0Preconditioner()
    except AttributeError:
        # If we are using a direct solver Pressure_Preconditioner_type will not
        # exist. This is not an error.
        press_precond = None

    # Pressure solver
    if ProjectParameters.Pressure_Linear_Solver == "Conjugate gradient":
        if ProjectParameters.Pressure_Preconditioner_type == 'None':
            pressure_linear_solver = CGSolver(ProjectParameters.Pressure_Iterative_Tolerance, ProjectParameters.Pressure_Solver_Max_Iteration)
        else:
            pressure_linear_solver = CGSolver(ProjectParameters.Pressure_Iterative_Tolerance, ProjectParameters.Pressure_Solver_Max_Iteration, press_precond)
    elif ProjectParameters.Pressure_Linear_Solver == "BiConjugate gradient stabilized":
        if ProjectParameters.Pressure_Preconditioner_type == 'None':
            pressure_linear_solver = BICGSTABSolver(ProjectParameters.Pressure_Iterative_Tolerance, ProjectParameters.Pressure_Solver_Max_Iteration)
        else:
            pressure_linear_solver = BICGSTABSolver(ProjectParameters.Pressure_Iterative_Tolerance, ProjectParameters.Pressure_Solver_Max_Iteration, press_precond)
    elif ProjectParameters.Pressure_Linear_Solver == "GMRES":
        if ProjectParameters.Pressure_Preconditioner_type == 'None':
            pressure_linear_solver = GMRESSolver(ProjectParameters.Pressure_Iterative_Tolerance, ProjectParameters.Pressure_Solver_Max_Iteration)
        else:
            pressure_linear_solver = GMRESSolver(ProjectParameters.Pressure_Iterative_Tolerance, ProjectParameters.Pressure_Solver_Max_Iteration, press_precond)
    elif ProjectParameters.Pressure_Linear_Solver == "Skyline LU factorization":
        pressure_linear_solver = SkylineLUFactorizationSolver()
    elif ProjectParameters.Pressure_Linear_Solver == "Super LU":
        pressure_linear_solver = SuperLUSolver()
    elif ProjectParameters.Pressure_Linear_Solver == "Parallel MKL Pardiso":
        from KratosMultiphysics.MKLSolversApplication import MKLPardisoSolver
        pressure_linear_solver = MKLPardisoSolver()

# copy Y_WALL
for node in fluid_model_part.Nodes:
    y = node.GetSolutionStepValue(Y_WALL, 0)
    node.SetValue(Y_WALL, y)

dynamic_tau = ProjectParameters.use_dt_in_stabilization
oss_switch = ProjectParameters.use_orthogonal_subscales
# creating the solvers
# fluid solver
if(ProjectParameters.SolverType == "FractionalStep"):
    fluid_solver = solver.IncompressibleFluidSolver(fluid_model_part, domain_size)
    fluid_solver.max_val_its = ProjectParameters.max_vel_its
    fluid_solver.max_press_its = ProjectParameters.max_press_its
    fluid_solver.predictor_corrector = ProjectParameters.predictor_corrector
    fluid_solver.vel_toll = ProjectParameters.velocity_relative_tolerance
    fluid_solver.press_toll = ProjectParameters.pressure_relative_tolerance
    fluid_solver.dynamic_tau = float(dynamic_tau)
    fluid_solver.compute_reactions = ProjectParameters.Calculate_reactions

fluid_solver.Initialize()
print("fluid solver created")

mesh_name = 0.0
gid_io.InitializeMesh(mesh_name)
gid_io.WriteMesh(fluid_model_part.GetMesh())
gid_io.FinalizeMesh()

# Stepping and time settings
full_Dt = ProjectParameters.Dt
initial_Dt = 0.001 * full_Dt  # 0.05 #0.01
Nsteps = ProjectParameters.nsteps

time = ProjectParameters.Start_time

#
import empire_wrapper

print('Create object wrapper ...')
wrapper = empire_wrapper.EmpireWrapper(fluid_model_part)

# loop information
it = 0
maxit = 1

for i in range(1, Nsteps + 1):
    print('################## Current time step: {}'.format(i))

    # while it < maxit:
    # print '--------------Iteration: {}'.format(it)

    if(i < 3):
        deltaT = initial_Dt

    else:
        deltaT = full_Dt

    time = time + deltaT
    fluid_model_part.CloneTimeStep(time)

    print("STEP = ", i)
    print("TIME = ", time)

    if(i >= 3):
        print('Receiving ...')
        wrapper.recvDisplacement()

        print('Solving ...')
        fluid_solver.Solve()

        print('Sending mesh ...')
        # wrapper.sendMesh()

        print('Sending forces ...')
        wrapper.sendForces()

        isConvergent = wrapper.recvConvergenceSignal()
        print('isConvergent: {}'.format(isConvergent))
        if isConvergent != 0:
            break

        gid_io.WriteNodalResults(PRESSURE, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(VELOCITY, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(MESH_VELOCITY, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(DISPLACEMENT, fluid_model_part.Nodes, time, 0)

            # it = it + 1

        # it = 1

gid_io.FinalizeResults()

#
# EMPIRE disconnect
# libempire_api.EMPIRE_API_Disconnect()
# exit();
