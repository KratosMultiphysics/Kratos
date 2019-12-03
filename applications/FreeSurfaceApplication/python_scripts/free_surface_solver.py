from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import problem_settings

# setting the domain size for the problem to be solved
domain_size = problem_settings.domain_size

# importing Kratos main library
import KratosMultiphysics
import KratosMultiphysics.FreeSurfaceApplication as KratosFreeSurface

# defining a model part for the fluid and one for the structure
FreeSurfaceModel = KratosMultiphysics.Model()
fluid_model_part = FreeSurfaceModel.CreateModelPart("FluidPart")

# importing the solvers needed
from KratosMultiphysics.FreeSurfaceApplication.edgebased_levelset_solver import EdgeBasedLevelSetSolver
EdgeBasedLevelSetSolver.AddVariables(fluid_model_part)

# introducing input file name
input_file_name = problem_settings.problem_name

# reading the fluid part
gid_mode = KratosMultiphysics.GiDPostMode.GiD_PostBinary
multifile = KratosMultiphysics.MultiFileFlag.MultipleFiles
deformed_mesh_flag = KratosMultiphysics.WriteDeformedMeshFlag.WriteUndeformed
write_conditions = KratosMultiphysics.WriteConditionsFlag.WriteConditions

# selecting output format
if problem_settings.print_layers:
    gid_io = KratosFreeSurface.EdgebasedGidIO(input_file_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)
else:
    gid_io = KratosMultiphysics.GidIO(input_file_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)

model_part_io_fluid = KratosMultiphysics.ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(2)

# adding dofs
EdgeBasedLevelSetSolver.AddDofs(fluid_model_part)

# we assume here that all of the internal nodes are marked with a negative distance
# set the distance of all of the internal nodes to a small value
small_value = 0.0001
n_active = 0
for node in fluid_model_part.Nodes:
    dist = node.GetSolutionStepValue( KratosMultiphysics.DISTANCE)
    if(dist < 0.0):
        n_active = n_active + 1
        node.SetSolutionStepValue( KratosMultiphysics.DISTANCE, 0, -small_value)
    else:
        node.SetSolutionStepValue( KratosMultiphysics.DISTANCE, 0, small_value)

if(n_active == 0):
    raise "ERROR. At least one node has to be initialized with a distance lesser than 0"

# make sure that the porosity is not zero on any node (set by default to fluid only)
for node in fluid_model_part.Nodes:
    if(node.GetSolutionStepValue( KratosMultiphysics.POROSITY) == 0.0):
        node.SetSolutionStepValue( KratosMultiphysics.POROSITY, 0, 1.0)
    if(node.GetSolutionStepValue( KratosMultiphysics.DIAMETER) == 0.0):
        node.SetSolutionStepValue( KratosMultiphysics.DIAMETER, 0, 1.0)

# constructing the solver
body_force = KratosMultiphysics.Vector(3)
body_force[0] = problem_settings.body_force_x
body_force[1] = problem_settings.body_force_y
body_force[2] = problem_settings.body_force_z
if(body_force[0] == 0.0 and body_force[1] == 0.0 and body_force[2] == 0.0):
    raise "ERROR. Body Force cannot be a ZERO VECTOR"

viscosity = problem_settings.viscosity
density = problem_settings.density
fluid_solver = EdgeBasedLevelSetSolver(fluid_model_part, domain_size, body_force, viscosity, density)
fluid_solver.redistance_frequency = problem_settings.redistance_frequency
fluid_solver.extrapolation_layers = int(problem_settings.extrapolation_layers)
fluid_solver.stabdt_pressure_factor = problem_settings.stabdt_pressure_factor
fluid_solver.stabdt_convection_factor = problem_settings.stabdt_convection_factor
fluid_solver.use_mass_correction = problem_settings.use_mass_correction
fluid_solver.tau2_factor = problem_settings.tau2_factor
fluid_solver.edge_detection_angle = problem_settings.edge_detection_angle
fluid_solver.assume_constant_pressure = problem_settings.assume_constant_pressure
fluid_solver.compute_porous_resistance_law = int(problem_settings.compute_porous_resistance_law)  # 0 = None; 1 = Ergun; 2 = Custom;
# print "compute_porous_resistance_law   ", fluid_solver.compute_porous_resistance_law
# using MKLPardisosolver ----> it has to be compiled in kratos!!
# fluid_solver.pressure_linear_solver = MKLPardisoSolver()

fluid_solver.Initialize()

if(problem_settings.wall_law_y > 1e-10):
    fluid_solver.fluid_solver.ActivateWallResistance(problem_settings.wall_law_y)

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

# mesh to be printed
if problem_settings.single_output_file:
    mesh_name = 0.0
    gid_io.InitializeMesh(mesh_name)
    gid_io.WriteMesh(fluid_model_part.GetMesh())
    gid_io.FinalizeMesh()
    gid_io.Flush()

    gid_io.InitializeResults(mesh_name, (fluid_model_part).GetMesh())

max_safety_factor = safety_factor

time = 0.0
step = 0
next_output_time = output_dt
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
            fluid_solver.fluid_solver.ReduceTimeStep(fluid_model_part, time)  # this is to set the database to the value at the beginning of the step

            safety_factor *= problem_settings.reduction_on_failure
            reduced_dt = fluid_solver.EstimateTimeStep(safety_factor, max_Dt)

            print("time before reduction= ", time)
            time = time - Dt + reduced_dt
            print("reduced time = ", time)
            print("Dt = ", Dt)
            print("reduced_dt = ", reduced_dt)

            fluid_solver.fluid_solver.ReduceTimeStep(fluid_model_part, time)  # this is to set the database to the value at the beginning of the step

            fluid_solver.Solve()

    if(time >= next_output_time):
        if not problem_settings.single_output_file:
            # writing mesh
            gid_io.InitializeMesh(time)
            gid_io.WriteMesh((fluid_model_part).GetMesh())
            gid_io.FinalizeMesh()
            gid_io.InitializeResults(time, (fluid_model_part).GetMesh())

        gid_io.WriteNodalResults( KratosMultiphysics.PRESSURE, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults( KratosMultiphysics.POROSITY, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults( KratosMultiphysics.VELOCITY, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults( KratosMultiphysics.DISTANCE, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults( KratosMultiphysics.PRESS_PROJ, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults( KratosMultiphysics.LIN_DARCY_COEF, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(KM.NONLIN_DARCY_COEF, fluid_model_part.Nodes, time, 0)
        gid_io.Flush()

        if not problem_settings.single_output_file:
            gid_io.FinalizeResults()

        next_output_time = time + output_dt

        out = 0

    out = out + 1
    step = step + 1

if problem_settings.single_output_file:
    gid_io.FinalizeResults()
