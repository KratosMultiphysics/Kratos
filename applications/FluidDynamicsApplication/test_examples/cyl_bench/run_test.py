from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#
#
# import the configuration data as read from the GiD
import ProjectParameters
import define_output


# importing benchmarking
kratos_benchmarking_path = '../../../../benchmarking'  # kratos_root/benchmarking
import sys
sys.path.append(kratos_benchmarking_path)
import benchmarking


def PrintResults(model_part):
    print("Writing results. Please run Gid for viewing results of analysis.")
    for variable_name in ProjectParameters.nodal_results:
        gid_io.WriteNodalResults(variables_dictionary[variable_name], model_part.Nodes, time, 0)
    for variable_name in ProjectParameters.gauss_points_results:
        gid_io.PrintOnGaussPoints(variables_dictionary[variable_name], model_part, time)

    if(ProjectParameters.TurbulenceModel == "Spalart-Allmaras"):
        gid_io.WriteNodalResults(VISCOSITY, model_part.Nodes, time, 0)
        # gid_io.WriteNodalResults(MOLECULAR_VISCOSITY,model_part.Nodes,time,0)
        # gid_io.WriteNodalResults(TURBULENT_VISCOSITY,model_part.Nodes,time,0)
        # gid_io.WriteNodalResults(DISTANCE,model_part.Nodes,time,0)


#
#
# setting the domain size for the problem to be solved
domain_size = ProjectParameters.domain_size

#
#
import sys
sys.path.append(ProjectParameters.kratos_path)
from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *

# defining variables to be used

variables_dictionary = {"PRESSURE": PRESSURE,
                        "VELOCITY": VELOCITY,
                        "REACTION": REACTION,
                        "DISTANCE": DISTANCE, }

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
    import vms_fractional_step_solver as solver
    solver.AddVariables(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian"):
    import vms_monolithic_solver as solver
    solver.AddVariables(fluid_model_part)
else:
    raise NameError("solver type not supported: options are FractionalStep  - monolithic_solver_eulerian")

if(ProjectParameters.TurbulenceModel == "Spalart-Allmaras"):
    fluid_model_part.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY)
    fluid_model_part.AddNodalSolutionStepVariable(MOLECULAR_VISCOSITY)
    fluid_model_part.AddNodalSolutionStepVariable(TEMP_CONV_PROJ)

    fluid_model_part.AddNodalSolutionStepVariable(DISTANCE)
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

if(ProjectParameters.VolumeOutput):
    if ProjectParameters.GiDWriteConditionsFlag:
        write_conditions = WriteConditionsFlag.WriteConditions
    else:
        write_conditions = WriteConditionsFlag.WriteElementsOnly
else:
    write_conditions = WriteConditionsFlag.WriteConditions

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
fluid_model_part.SetBufferSize(3)

# adding dofs
if(SolverType == "FractionalStep"):
    solver.AddDofs(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian"):
    solver.AddDofs(fluid_model_part)

if(ProjectParameters.TurbulenceModel == "Spalart-Allmaras"):
    for node in fluid_model_part.Nodes:
        node.AddDof(TURBULENT_VISCOSITY)

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

elif "monolithic_solver_eulerian":  # single coupled solver
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
    elif ProjectParameters.Monolithic_Linear_Solver == "MixedUP":
        velocity_linear_solver = SuperLUIterativeSolver()
        pressure_linear_solver = ScalingSolver(BICGSTABSolver(1e-3, 500), True)
        m = 5
        max_it = m
        monolithic_linear_solver = MixedUPLinearSolver(velocity_linear_solver, pressure_linear_solver, 1e-6, max_it, m)
    elif ProjectParameters.Monolithic_Linear_Solver == "SuperLUIterativeSolver":
        monolithic_linear_solver = ScalingSolver(SuperLUIterativeSolver(), True)

# copy Y_WALL
for node in fluid_model_part.Nodes:
    y = node.GetSolutionStepValue(Y_WALL, 0)
    node.SetValue(Y_WALL, y)

dynamic_tau = ProjectParameters.use_dt_in_stabilization
oss_switch = ProjectParameters.use_orthogonal_subscales
# creating the solvers
# fluid solver
if(SolverType == "FractionalStep"):
    fluid_solver = solver.IncompressibleFluidSolver(fluid_model_part, domain_size)
    fluid_solver.max_val_its = ProjectParameters.max_vel_its
    fluid_solver.max_press_its = ProjectParameters.max_press_its
    fluid_solver.predictor_corrector = ProjectParameters.predictor_corrector
    fluid_solver.vel_toll = ProjectParameters.velocity_relative_tolerance
    fluid_solver.press_toll = ProjectParameters.pressure_relative_tolerance
    fluid_solver.dynamic_tau = float(dynamic_tau)
    fluid_solver.compute_reactions = ProjectParameters.Calculate_reactions
elif(SolverType == "monolithic_solver_eulerian"):
    fluid_solver = solver.MonolithicSolver(fluid_model_part, domain_size)
    fluid_solver.oss_switch = int(oss_switch)
    fluid_solver.dynamic_tau = float(dynamic_tau)
    fluid_solver.rel_vel_tol = ProjectParameters.velocity_relative_tolerance
    fluid_solver.abs_vel_tol = ProjectParameters.velocity_absolute_tolerance
    fluid_solver.rel_pres_tol = ProjectParameters.pressure_relative_tolerance
    fluid_solver.abs_pres_tol = ProjectParameters.pressure_absolute_tolerance
    fluid_solver.max_iter = ProjectParameters.max_iterations
    fluid_solver.compute_reactions = ProjectParameters.Calculate_reactions
    fluid_solver.linear_solver = monolithic_linear_solver


# activate turbulence model
if(ProjectParameters.TurbulenceModel == "Smagorinsky-Lilly"):
    fluid_solver.activate_smagorinsky(ProjectParameters.SmagorinskyConstant)
elif(ProjectParameters.TurbulenceModel == "Spalart-Allmaras"):
    # apply the initial turbulent viscosity on all of the nodes
    turb_visc = ProjectParameters.TurbulentViscosity
    for node in fluid_model_part.Nodes:
        node.SetSolutionStepValue(TURBULENT_VISCOSITY, 0, turb_visc);
        visc = node.GetSolutionStepValue(VISCOSITY)
        node.SetSolutionStepValue(MOLECULAR_VISCOSITY, 0, visc);
        if(node.IsFixed(VELOCITY_X)):
            node.Fix(TURBULENT_VISCOSITY)

    # select nodes on the wall
    wall_nodes = []
    for i in ProjectParameters.SA_wall_group_ids:
        nodes = fluid_model_part.GetNodes(i)  # get the nodes of the wall for SA.
        for node in nodes:
            wall_nodes.append(node)
            node.SetSolutionStepValue(TURBULENT_VISCOSITY, 0, 0.0);
            node.Fix(TURBULENT_VISCOSITY)

    DES = False
    fluid_solver.ActivateSpalartAllmaras(wall_nodes, DES)


fluid_solver.Initialize()


print("fluid solver created")


cut_model_part = model.CreateModelPart("CutPart");
if(ProjectParameters.VolumeOutput):
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
else:
    # generate the cuts
    Cut_App = Cutting_Application()
    Cut_App.FindSmallestEdge(fluid_model_part)
    Cut_App.AddVariablesToCutModelPart(fluid_model_part,cut_model_part)

    cut_number = 1

    cut_list = define_output.DefineCutPlanes()
    print(cut_list)

    for item in cut_list:
        print(item)
        n = Vector(item[0])
        p = Vector(item[1])
        Cut_App.GenerateCut(fluid_model_part, cut_model_part, n, p, cut_number, 0.01)
        cut_number = cut_number + 1
        print("generated cut number =", cut_number)

    Cut_App.AddSkinConditions(fluid_model_part, cut_model_part, cut_number)
    cut_number += 1

    # mesh to be printed (single mesh case)
    if ProjectParameters.GiDMultiFileFlag == "Single":
        mesh_name = 0.0
        gid_io.InitializeMesh(mesh_name)
        gid_io.WriteMesh(cut_model_part.GetMesh())
        gid_io.FinalizeMesh()

        gid_io.InitializeResults(mesh_name, (cut_model_part).GetMesh())
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

# 33
# 33
# define the drag computation list
drag_list = define_output.DefineDragList()
drag_file_output_list = []
for it in drag_list:
    f = open(it[1], 'w')
    drag_file_output_list.append(f)
    tmp = "#Drag for group " + it[1] + "\n"
    f.write(tmp)
    tmp = "time RX RY RZ"
    f.write(tmp)
    f.flush()

print(drag_file_output_list)


def PrintDrag(drag_list, drag_file_output_list, fluid_model_part, time):
    i = 0
    for it in drag_list:
        print(it[0])
        nodes = fluid_model_part.GetNodes(it[0])
        drag = Vector(3);
        drag[0] = 0.0
        drag[1] = 0.0
        drag[2] = 0.0
        for node in nodes:
            reaction = node.GetSolutionStepValue(REACTION, 0)
            drag[0] += reaction[0]
            drag[1] += reaction[1]
            drag[2] += reaction[2]

        if (benchmarking.InBenchmarkingMode()):
            benchmarking.Output(drag[0], "drag X", 1e-9, 1e-3)
            benchmarking.Output(drag[1], "drag Y", 1e-9, 1e-3)
            benchmarking.Output(drag[2], "drag Z", 1e-9, 1e-3)

        output = str(time) + " " + str(drag[0]) + " " + str(drag[1]) + " " + str(drag[2]) + "\n"
        # print drag_file_output_list[i]
        # print output
        drag_file_output_list[i].write(output)
        drag_file_output_list[i].flush()
        i = i + 1


# 33
# preparing output of point graphs
import point_graph_printer

output_nodes_list = define_output.DefineOutputPoints()
graph_printer = point_graph_printer.PrintGraphPrinter(output_nodes_list, fluid_model_part, variables_dictionary, domain_size)


# Stepping and time settings
Dt = ProjectParameters.Dt
full_Dt = Dt
initial_Dt = 0.001 * full_Dt  # 0.05 #0.01
Nsteps = ProjectParameters.nsteps
final_time = ProjectParameters.max_time
output_time = ProjectParameters.output_time

time = ProjectParameters.Start_time
out = 0
step = 0

while(time <= final_time):

    if(step < 3):
        Dt = initial_Dt
    else:
        Dt = full_Dt

    time = time + Dt
    step = step + 1
    fluid_model_part.CloneTimeStep(time)

    print("STEP = ", step)
    print("TIME = ", time)

    if(step >= 3):
        fluid_solver.Solve()

        if(step < 4):
            for k in range(0, ProjectParameters.divergence_cleareance_step):
                print("DOING DIVERGENCE CLEAREANCE")
                buffer_size = fluid_model_part.GetBufferSize()
                for i in range(0, buffer_size):
                    for node in fluid_model_part.Nodes:
                        vel = node.GetSolutionStepValue(VELOCITY)
                        node.SetSolutionStepValue(VELOCITY, i, vel)
                        node.SetSolutionStepValue(PRESSURE, i, 0.0)
                    if(SolverType == "monolithic_solver_eulerian"):
                        zero_vector = Vector(3)
                        zero_vector[0] = 0.0
                        zero_vector[1] = 0.0
                        zero_vector[2] = 0.0
                        for node in fluid_model_part.Nodes:
                            node.SetSolutionStepValue(ACCELERATION, i, zero_vector)
                    if(ProjectParameters.TurbulenceModel == "Spalart-Allmaras"):
                        for node in fluid_model_part.Nodes:
                            visc = node.GetSolutionStepValue(VISCOSITY)
                            node.SetSolutionStepValue(VISCOSITY, i, visc)

                fluid_solver.Solve()

        graph_printer.PrintGraphs(time)
        PrintDrag(drag_list, drag_file_output_list, fluid_model_part, time)

    if(output_time <= out):
        if(ProjectParameters.VolumeOutput):
            if Multifile:
                gid_io.InitializeMesh(time)
                gid_io.WriteMesh(fluid_model_part.GetMesh())
                gid_io.FinalizeMesh()

                gid_io.InitializeResults(time, (fluid_model_part).GetMesh())

            PrintResults(fluid_model_part)
            out = 0
        else:
            cut_model_part.CloneTimeStep(time)
            Cut_App.UpdateCutData(cut_model_part, fluid_model_part)

            if Multifile:
                gid_io.InitializeMesh(time)
                gid_io.WriteMesh(cut_model_part.GetMesh())
                gid_io.FinalizeMesh()

                gid_io.InitializeResults(time, (cut_model_part).GetMesh())

            PrintResults(cut_model_part)

            out = 0

        if Multifile:
            gid_io.FinalizeResults()
            if ProjectParameters.GiDPostMode == "Binary":
                f.write(ProjectParameters.problem_name + '_' + str(time) + '.post.bin\n')
            elif ProjectParameters.GiDPostMode == "Ascii":
                f.write(ProjectParameters.problem_name + '_' + str(time) + '.post.msh\n')

    out = out + Dt

if Multifile:
    f.close()
else:
    gid_io.FinalizeResults()

for i in drag_file_output_list:
    i.close();
