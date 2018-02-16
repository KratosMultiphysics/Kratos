from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import problem_settings

#
#
# setting the domain size for the problem to be solved
domain_size = problem_settings.domain_size

#
#
# ATTENTION: here the order is important

# including kratos path
import sys
sys.path.append(problem_settings.kratos_path)

# importing Kratos main library
from KratosMultiphysics import *

# from now on the order is not anymore crucial
#
#
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.MeshingApplication import *

# defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart")
cut_model_part = ModelPart("CutPart")

#


# importing the solvers needed
import edgebased_levelset_solver
edgebased_levelset_solver.AddVariables(fluid_model_part)

# introducing input file name
input_file_name = problem_settings.problem_name

# reading the fluid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteConditions

# selecting output format
gid_io = GidIO(
    input_file_name,
    gid_mode,
    multifile,
    deformed_mesh_flag,
    write_conditions)

model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(2)

# adding dofs
edgebased_levelset_solver.AddDofs(fluid_model_part)

# we assume here that all of the internal nodes are marked with a negative distance
# set the distance of all of the internal nodes to a small value
small_value = 0.0001
n_active = 0
for node in fluid_model_part.Nodes:
    dist = node.GetSolutionStepValue(DISTANCE)
    if(dist < 0.0):
        n_active = n_active + 1
        node.SetSolutionStepValue(DISTANCE, 0, -small_value)
    else:
        node.SetSolutionStepValue(DISTANCE, 0, small_value)

if(n_active == 0):
    raise "ERROR. At least one node has to be initialized with a distance lesser than 0"

# make sure that the porosity is not zero on any node (set by default to
# fluid only)
for node in fluid_model_part.Nodes:
    if(node.GetSolutionStepValue(POROSITY) == 0.0):
        node.SetSolutionStepValue(POROSITY, 0, 1.0)
    if(node.GetSolutionStepValue(DIAMETER) == 0.0):
        node.SetSolutionStepValue(DIAMETER, 0, 1.0)

# constructing the solver
body_force = Vector(3)
body_force[0] = problem_settings.body_force_x
body_force[1] = problem_settings.body_force_y
body_force[2] = problem_settings.body_force_z
if(body_force[0] == 0.0 and body_force[1] == 0.0 and body_force[2] == 0.0):
    raise "ERROR. Body Force cannot be a ZERO VECTOR"

viscosity = problem_settings.viscosity
density = problem_settings.density
fluid_solver = edgebased_levelset_solver.EdgeBasedLevelSetSolver(
    fluid_model_part, domain_size, body_force, viscosity, density)
fluid_solver.redistance_frequency = problem_settings.redistance_frequency
fluid_solver.extrapolation_layers = int(problem_settings.extrapolation_layers)
fluid_solver.stabdt_pressure_factor = problem_settings.stabdt_pressure_factor
fluid_solver.stabdt_convection_factor = problem_settings.stabdt_convection_factor
fluid_solver.use_mass_correction = problem_settings.use_mass_correction
fluid_solver.tau2_factor = problem_settings.tau2_factor
fluid_solver.edge_detection_angle = problem_settings.edge_detection_angle
fluid_solver.assume_constant_pressure = problem_settings.assume_constant_pressure
# 0 = None; 1 = Ergun; 2 = Custom;
fluid_solver.compute_porous_resistance_law = int(
    problem_settings.compute_porous_resistance_law)

pressure_fixed = problem_settings.pressure_fixed
fix_location = problem_settings.fix_location
Z_coord_free_surface = problem_settings.Z_coord_free_surface
Z_coord_bottom = problem_settings.Z_coord_bottom
X_coord_I_O = problem_settings.X_coord_I_O
free_surface_Z_coord = problem_settings.free_surface_Z_coord
X1 = problem_settings.X1
Y1 = problem_settings.Y1
X2 = problem_settings.X2
Y2 = problem_settings.Y2
X3 = problem_settings.X3
Y3 = problem_settings.Y3
X4 = problem_settings.X4
Y4 = problem_settings.Y4
X5 = problem_settings.X5
Y5 = problem_settings.Y5
Xtol = problem_settings.Xtol
Ytol = problem_settings.Ytol
# print "compute_porous_resistance_law   ", fluid_solver.compute_porous_resistance_law

if(pressure_fixed == "ON"):
    node_list = []
    for node in fluid_model_part.Nodes:
        if(fix_location == "outlet"):
            if(node.X > X_coord_I_O - 0.001):
                if(node.Z < Z_coord_free_surface + 0.001):
                    if(node.Z > Z_coord_bottom - 0.001):
                        node_list.append(node)
        if(fix_location == "inlet"):
            if(node.X < X_coord_I_O + 0.001):
                if(node.Z < Z_coord_free_surface + 0.001):
                    if(node.Z > Z_coord_bottom - 0.001):
                        node_list.append(node)

    import math

    for node in node_list:
        Zcoord = node.Z
        p = (Z_coord_free_surface - Zcoord) * (-body_force[2]) * density
        node.SetSolutionStepValue(PRESSURE, 0, p)
        node.Fix(PRESSURE)

fluid_solver.Initialize()

if(problem_settings.wall_law_y > 1e-10):
    fluid_solver.fluid_solver.ActivateWallResistance(
        problem_settings.wall_law_y)

#


print("fluid solver created")

# settings to be changed
max_Dt = problem_settings.max_time_step
initial_Dt = 0.001 * max_Dt
final_time = problem_settings.max_time
output_dt = problem_settings.output_dt
safety_factor = problem_settings.safety_factor

number_of_inital_steps = problem_settings.number_of_inital_steps
initial_time_step = problem_settings.initial_time_step
out = 0

original_max_dt = max_Dt

max_safety_factor = safety_factor

time = 0.0
step = 0
next_output_time = output_dt


if(free_surface_Z_coord == 1):
    number_points_x = 1
    number_points_y = 1
    X_position = [X1]
    Y_position = [Y1]
    position_tolerance_X = Xtol
    position_tolerance_Y = Ytol

if(free_surface_Z_coord == 2):
    number_points_x = 2
    number_points_y = 2
    X_position = [X1, X2]
    Y_position = [Y1, Y2]
    position_tolerance_X = Xtol
    position_tolerance_Y = Ytol

if(free_surface_Z_coord == 3):
    number_points_x = 3
    number_points_y = 3
    X_position = [X1, X2, X3]
    Y_position = [Y1, Y2, Y3]
    position_tolerance_X = Xtol
    position_tolerance_Y = Ytol

if(free_surface_Z_coord == 4):
    number_points_x = 4
    number_points_y = 4
    X_position = [X1, X2, X3, X4]
    Y_position = [Y1, Y2, Y3, Y4]
    position_tolerance_X = Xtol
    position_tolerance_Y = Ytol

if(free_surface_Z_coord == 5):
    number_points_x = 5
    number_points_y = 5
    X_position = [X1, X2, X3, X4, X5]
    Y_position = [Y1, Y2, Y3, Y4, Y5]
    position_tolerance_X = Xtol
    position_tolerance_Y = Ytol

if(free_surface_Z_coord > 0):
    f = open("altura_sup_lib.txt", "w")
    f.write('     Tiempo          ')
    i = 0
    for j in range(0, (number_points_y)):
        for i in range(0, (number_points_x)):
            if(i == j):
                if (i + j < (number_points_x + number_points_y - 2)):
                    f.write(
                        'X = ' + str(X_position[i]) + ', Y = ' + str(Y_position[j]) + '    ')
                if (i + j == (number_points_x + number_points_y - 2)):
                    f.write(
                        'X = ' + str(X_position[i]) + ', Y = ' + str(Y_position[j]) + '\n')
            i = i + 1
        j = j + 1
    f.close()

while(time < final_time):

    if(step < number_of_inital_steps):
        max_Dt = initial_time_step
    else:
        max_Dt = original_max_dt
        # progressively increment the safety factor
        # in the steps that follow a reduction of it
        safety_factor = safety_factor * 1.2
        if(safety_factor > max_safety_factor):
            safety_factor = max_safety_factor

    Dt = fluid_solver.EstimateTimeStep(safety_factor, max_Dt)

    time = time + Dt
    fluid_model_part.CloneTimeStep(time)
    cut_model_part.CloneTimeStep(time)

    print("******** CURRENT TIME = ", time)

    if(step >= 3):
        fluid_solver.Solve()

        check_dt = fluid_solver.EstimateTimeStep(0.95, max_Dt)

        if(check_dt < Dt):
            print("***********************************************************")
            print("***********************************************************")
            print("***********************************************************")
            print("            *** REDUCING THE TIME STEP ***")
            print("***********************************************************")
            print("***********************************************************")
            print("***********************************************************")

            # we found a velocity too large! we need to reduce the time step
            # this is to set the database to the value at the beginning of the
            # step
            fluid_solver.fluid_solver.ReduceTimeStep(fluid_model_part, time)

            safety_factor *= problem_settings.reduction_on_failure
            reduced_dt = fluid_solver.EstimateTimeStep(safety_factor, max_Dt)

            print("time before reduction= ", time)
            time = time - Dt + reduced_dt
            print("reduced time = ", time)
            print("Dt = ", Dt)
            print("reduced_dt = ", reduced_dt)

            # this is to set the database to the value at the beginning of the
            # step
            fluid_solver.fluid_solver.ReduceTimeStep(fluid_model_part, time)

            fluid_solver.Solve()

    if(time >= next_output_time):

        Cut_App = Cutting_Isosurface_Application()
        Cut_App.DeleteCutData(cut_model_part)
        variable = DISTANCE
        isovalue = 0
        tolerance = 0.000000001
        Cut_App.GenerateScalarVarCut(
            fluid_model_part,
            cut_model_part,
            variable,
            isovalue,
            1,
            tolerance)

        if(free_surface_Z_coord > 0):

            f = open("altura_sup_lib.txt", "a")
            f.write(str(time) + '      ')

            Cut_App.UpdateCutData(cut_model_part, fluid_model_part)

            i = 1
            for j in range(0, (number_points_y)):
                for i in range(0, (number_points_x)):
                    n = 0
                    h = 0
                    if(i == j):
                        for node in cut_model_part.Nodes:
                            if (node.X > (X_position[i] - position_tolerance_X) and node.X < (X_position[i] + position_tolerance_X)):
                                if (node.Y > (Y_position[j] - position_tolerance_Y) and node.Y < (Y_position[j] + position_tolerance_Y)):
                                    n = n + 1
                                    h = h + node.Z
                        if (i + j < (number_points_x + number_points_y - 2)):
                            if (n > 0):
                                h_tot = h / n
                                print("Altura de la superficie libre en X = ", X_position[i], " , Y = ", Y_position[j], " , ", h_tot)
                                f.write(str(h_tot) + '      ')
                            else:
                                f.write('No data available      ')
                        if (i + j == (number_points_x + number_points_y - 2)):
                            if (n > 0):
                                h_tot = h / n
                                print("Altura de la superficie libre en X = ", X_position[i], " , Y = ", Y_position[j], " , ", h_tot)
                                f.write(str(h_tot) + '      ' + '\n')
                            else:
                                f.write('No data available      ' + '\n')
                    i = i + 1
                j = j + 1
            f.close()

        Cut_App.AddModelPartElements(fluid_model_part, cut_model_part, 2)

        # meh to be printed
        gid_io.InitializeMesh(time)
        gid_io.WriteMesh((cut_model_part).GetMesh())
        gid_io.FinalizeMesh()
        gid_io.InitializeResults(time, (cut_model_part).GetMesh())
        Cut_App.UpdateCutData(cut_model_part, fluid_model_part)

        gid_io.WriteNodalResults(PRESSURE, cut_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(POROSITY, cut_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(VELOCITY, cut_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(DISTANCE, cut_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(PRESS_PROJ, cut_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(LIN_DARCY_COEF, cut_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(
            NONLIN_DARCY_COEF,
            cut_model_part.Nodes,
            time,
            0)
        gid_io.Flush()

        gid_io.FinalizeResults()

        next_output_time = time + output_dt

        out = 0

    out = out + 1
    step = step + 1

gid_io.FinalizeResults()
