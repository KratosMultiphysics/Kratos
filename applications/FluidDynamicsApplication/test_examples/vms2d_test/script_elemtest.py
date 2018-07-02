from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#
#
# import the configuration data as read from the GiD
import ProjectParameters


def PrintResults(model_part):
    print("Writing results. Please run Gid for viewing results of analysis.")
    for variable_name in ProjectParameters.nodal_results:
        gid_io.WriteNodalResults(variables_dictionary[variable_name], model_part.Nodes, time, 0)
    for variable_name in ProjectParameters.gauss_points_results:
        gid_io.PrintOnGaussPoints(variables_dictionary[variable_name], model_part, time)


def FindNode(node_list, x, y, z):
    for node in node_list:
        if ((node.X - x) ** 2 + (node.Y - y) ** 2 + (node.Z - z) ** 2 < 0.0000001):
            return node


def BenchmarkCheck(time, node1, node2):
    benchmarking.Output(time, "Time",1e-7)
    benchmarking.Output(node1.GetSolutionStepValue(PRESSURE), "Test node 1 pressure", None, 0.01)
    benchmarking.Output(node2.GetSolutionStepValue(VELOCITY_X), "Test node 2 velocity x", None, 0.01)

#
domain_size = ProjectParameters.domain_size

kratos_path = '../../../..'
kratos_benchmarking_path = '../../../../benchmarking'  # kratos_root/benchmarking
import sys
sys.path.append(kratos_path)
sys.path.append(kratos_benchmarking_path)

from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *

import benchmarking

# defining variables to be used

variables_dictionary = {"PRESSURE": PRESSURE,
                        "VELOCITY": VELOCITY,
                        "REACTION": REACTION,
                        "DISTANCE": DISTANCE,
                        "ADVPROJ": ADVPROJ,
                        "DIVPROJ": DIVPROJ,
                        "VORTICITY": VORTICITY, }

# defining a model part for the fluid
model = Model()
fluid_model_part = model.CreateModelPart("FluidPart")


if "REACTION" in ProjectParameters.nodal_results:
    fluid_model_part.AddNodalSolutionStepVariable(REACTION)
if "DISTANCE" in ProjectParameters.nodal_results:
    fluid_model_part.AddNodalSolutionStepVariable(DISTANCE)

#
# importing the solvers needed
SolverType = ProjectParameters.SolverType
if(SolverType == "FractionalStep"):
    import incompressible_fluid_solver
    incompressible_fluid_solver.AddVariables(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian"):
    import vms_monolithic_solver
    vms_monolithic_solver.AddVariables(fluid_model_part)
else:
    raise NameError("solver type not supported: options are FractionalStep - monolithic_solver_eulerian")

# introducing input file name
input_file_name = ProjectParameters.problem_name

# reading the fluid part

# initialize GiD  I/O
if ProjectParameters.GiDPostMode == "Binary":
    gid_mode = GiDPostMode.GiD_PostBinary
elif ProjectParameters.GiDPostMode == "Ascii":
    gid_mode = GiDPostMode.GiD_PostAscii
elif ProjectParameters.GiDPostMode == "AsciiZipped":
    gid_mode = GiDPostMode.GiD_PostAsciiZipped
else:
    print("Unknown GiD post mode, assuming Binary")
    gid_mode = GiDPostMode.GiD_PostBinary

if ProjectParameters.GiDWriteMeshFlag:
    deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
else:
    deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed

if ProjectParameters.GiDWriteConditionsFlag:
    write_conditions = WriteConditionsFlag.WriteConditions
else:
    write_conditions = WriteConditionsFlag.WriteElementsOnly

if ProjectParameters.GiDMultiFileFlag == "Single":
    multifile = MultiFileFlag.SingleFile
elif ProjectParameters.GiDMultiFileFlag == "Multiples":
    multifile = MultiFileFlag.MultipleFiles
else:
    print("Unknown GiD multiple file mode, assuming Single")
    multifile = MultiFileFlag.SingleFile

gid_io = GidIO(input_file_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
if SolverType == "FractionalStep":
    fluid_model_part.SetBufferSize(3)
else:
    fluid_model_part.SetBufferSize(2)

# adding dofs
if(SolverType == "FractionalStep"):
    incompressible_fluid_solver.AddDofs(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian"):
    vms_monolithic_solver.AddDofs(fluid_model_part)

# If Lalplacian form = 2, free all pressure Dofs
laplacian_form = ProjectParameters.laplacian_form
if(laplacian_form >= 2):
    for node in fluid_model_part.Nodes:
        node.Free(PRESSURE)

# Commented out for efficiency, problem type assigns the values automatically
# check to ensure that no node has zero density or pressure
# for node in fluid_model_part.Nodes:
# if(node.GetSolutionStepValue(DENSITY) == 0.0):
# print "node ",node.Id," has zero density!"
# raise 'node with zero density found'
# if(node.GetSolutionStepValue(VISCOSITY) == 0.0):
# print "node ",node.Id," has zero viscosity!"
# raise 'node with zero VISCOSITY found'

# decoupled solver schemes: need solver for pressure and solver fol velocity

if ProjectParameters.SolverType == "FractionalStep":
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
        
elif ProjectParameters.SolverType == "monolithic_solver_eulerian":  # single coupled solver
    # preconditioner
    try:
        if ProjectParameters.Monolithic_Preconditioner_type == 'Diagonal':
            precond = DiagonalPreconditioner()
        elif ProjectParameters.Monolithic_Preconditioner_type == 'ILU0':
            precond = ILU0Preconditioner()
    except AttributeError:
        # If we are using a direct solver Pressure_Preconditioner_type will not
        # exist. This is not an error.
        precond = None

    # solver
    if ProjectParameters.Monolithic_Linear_Solver == "Conjugate gradient":
        if ProjectParameters.Monolithic_Preconditioner_type == 'None':
            monolithic_linear_solver = CGSolver(ProjectParameters.Monolithic_Iterative_Tolerance, ProjectParameters.Monolithic_Solver_Max_Iteration)
        else:
            monolithic_linear_solver = CGSolver(ProjectParameters.Monolithic_Iterative_Tolerance, ProjectParameters.Monolithic_Solver_Max_Iteration, precond)
    elif ProjectParameters.Monolithic_Linear_Solver == "BiConjugate gradient stabilized":
        if ProjectParameters.Monolithic_Preconditioner_type == 'None':
            monolithic_linear_solver = BICGSTABSolver(ProjectParameters.Monolithic_Iterative_Tolerance, ProjectParameters.Monolithic_Solver_Max_Iteration)
        else:
            monolithic_linear_solver = BICGSTABSolver(ProjectParameters.Monolithic_Iterative_Tolerance, ProjectParameters.Monolithic_Solver_Max_Iteration, precond)
    elif ProjectParameters.Monolithic_Linear_Solver == "GMRES":
        if ProjectParameters.Monolithic_Preconditioner_type == 'None':
            monolithic_linear_solver = GMRESSolver(ProjectParameters.Monolithic_Iterative_Tolerance, ProjectParameters.Monolithic_Solver_Max_Iteration)
        else:
            monolithic_linear_solver = GMRESSolver(ProjectParameters.Monolithic_Iterative_Tolerance, ProjectParameters.Monolithic_Solver_Max_Iteration, precond)
    elif ProjectParameters.Monolithic_Linear_Solver == "Skyline LU factorization":
        monolithic_linear_solver = SkylineLUFactorizationSolver()
    elif ProjectParameters.Monolithic_Linear_Solver == "Super LU":
        monolithic_linear_solver = SuperLUSolver()

dynamic_tau = ProjectParameters.use_dt_in_stabilization
oss_switch = ProjectParameters.use_orthogonal_subscales
# creating the solvers
# fluid solver
if(SolverType == "FractionalStep"):
    fluid_solver = incompressible_fluid_solver.IncompressibleFluidSolver(fluid_model_part, domain_size)
    fluid_solver.max_val_its = ProjectParameters.max_vel_its
    fluid_solver.max_press_its = ProjectParameters.max_press_its
    fluid_solver.laplacian_form = laplacian_form
    # standard laplacian form
    fluid_solver.predictor_corrector = ProjectParameters.predictor_corrector
    fluid_solver.use_dt_in_stabilization = False
    fluid_solver.oss_switch = ProjectParameters.use_orthogonal_subscales
    fluid_solver.dynamic_tau = ProjectParameters.use_dt_in_stabilization
    fluid_solver.vel_toll = ProjectParameters.velocity_relative_tolerance
    fluid_solver.press_toll = ProjectParameters.pressure_relative_tolerance
    fluid_solver.CalculateReactions = ProjectParameters.Calculate_reactions
    # Solver definition
    fluid_solver.velocity_linear_solver = velocity_linear_solver
    fluid_solver.pressure_linear_solver = pressure_linear_solver
    fluid_solver.Initialize()
    
elif(SolverType == "monolithic_solver_eulerian"):
    fluid_solver = vms_monolithic_solver.MonolithicSolver(fluid_model_part, domain_size)
    fluid_solver.oss_switch = ProjectParameters.use_orthogonal_subscales
    fluid_solver.dynamic_tau = ProjectParameters.use_dt_in_stabilization
    fluid_solver.rel_vel_tol = ProjectParameters.velocity_relative_tolerance
    fluid_solver.abs_vel_tol = ProjectParameters.velocity_absolute_tolerance
    fluid_solver.rel_pres_tol = ProjectParameters.pressure_relative_tolerance
    fluid_solver.abs_pres_tol = ProjectParameters.pressure_absolute_tolerance
    fluid_solver.max_iter = ProjectParameters.max_iterations
    fluid_solver.CalculateReactions = ProjectParameters.Calculate_reactions
    # Solver definition
    fluid_solver.linear_solver = monolithic_linear_solver
    # fluid_solver.pressure_linear_solver = pressure_linear_solver
    fluid_solver.Initialize()
    
print("fluid solver created")

# mesh to be printed (single mesh case)
if ProjectParameters.GiDMultiFileFlag == "Single":
    mesh_name = 0.0
    gid_io.InitializeMesh(mesh_name)
    gid_io.WriteMesh(fluid_model_part.GetMesh())
    gid_io.FinalizeMesh()

    gid_io.InitializeResults(mesh_name, (fluid_model_part).GetMesh())
    Multifile = False

    # Write .post.list file (GiD postprocess list)
    f = open(ProjectParameters.problem_name + '.post.lst', 'w')
    f.write('Single\n')
    if ProjectParameters.GiDPostMode == "Binary":
        f.write(ProjectParameters.problem_name + '.post.bin\n')
    elif ProjectParameters.GiDPostMode == "Ascii":
        f.write(ProjectParameters.problem_name + '_0.post.msh\n')
    f.close()

else:  # ProjectParameters.GiDMultiFileFlag == "Multiples":
    Multifile = True

    # Initialize .post.list file (GiD postprocess list)
    f = open(ProjectParameters.problem_name + '.post.lst', 'w')
    f.write('Multiple\n')

# Benchmark nodes
node802 = FindNode(fluid_model_part.Nodes, 9.20572, 3.96404, 0.00000)  # Node 802
node1013 = FindNode(fluid_model_part.Nodes, 7.07261, 6.38485, 0.00000)  # Node 1013

# Stepping and time settings
Dt = ProjectParameters.Dt
full_Dt = Dt
initial_Dt = 0.01 * full_Dt  # 0.05 #0.01
Nsteps = ProjectParameters.nsteps
final_time = ProjectParameters.max_time
output_time = ProjectParameters.output_time

time = ProjectParameters.Start_time
out = 0
step = 0

while(time <= final_time):

    if(step < 10):
        Dt = initial_Dt
    else:
        Dt = full_Dt
 #   if step > 100:
 #       fluid_model_part.ProcessInfo.SetValue(OSS_SWITCH,1)

    time = time + Dt
    step = step + 1
    fluid_model_part.CloneTimeStep(time)

    print("STEP = ", step)
    print("TIME = ", time)

    if(step >= 3):
        fluid_solver.Solve()

    if(time >= 20):
        # Switching on OSS for test purposes
        fluid_model_part.ProcessInfo.SetValue(OSS_SWITCH, 1)

    if(output_time <= out):
        if Multifile:
            gid_io.InitializeMesh(time)
            gid_io.WriteMesh(fluid_model_part.GetMesh())
            gid_io.FinalizeMesh()

            gid_io.InitializeResults(time, (fluid_model_part).GetMesh())

        PrintResults(fluid_model_part)
        out = 0

        if Multifile:
            gid_io.FinalizeResults()
            if ProjectParameters.GiDPostMode == "Binary":
                f.write(ProjectParameters.problem_name + '_' + str(time) + '.post.bin\n')
            elif ProjectParameters.GiDPostMode == "Ascii":
                f.write(ProjectParameters.problem_name + '_' + str(time) + '.post.msh\n')

    out = out + Dt

    BenchmarkCheck(time, node802, node1013)

if Multifile:
    f.close()
else:
    gid_io.FinalizeResults()
