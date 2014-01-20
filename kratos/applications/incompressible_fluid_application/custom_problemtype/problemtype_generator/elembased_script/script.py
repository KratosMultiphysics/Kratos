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
from KratosMultiphysics.ConvectionDiffusionApplication import *

# defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart")

#


# importing the solvers needed
import level_set_elembased_fluid_solver
level_set_elembased_fluid_solver.AddVariables(fluid_model_part)

# introducing input file name
input_file_name = problem_settings.problem_name

# reading the fluid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly

# selecting output format
if(problem_settings.print_layers):
    gid_io = EdgebasedGidIO(input_file_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)
else:
    gid_io = GidIO(input_file_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)

model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(2)

# adding dofs
level_set_elembased_fluid_solver.AddDofs(fluid_model_part)

# we assume here that all of the internal nodes are marked with a negative distance
# set the distance of all of the internal nodes to a small value
small_value = 0.0001
n_active = 0
for node in fluid_model_part.Nodes:
    node.SetSolutionStepValue(VISCOSITY, 0, problem_settings.viscosity)
    node.SetSolutionStepValue(DENSITY, 0, problem_settings.density)
    node.SetSolutionStepValue(BODY_FORCE_X, 0, problem_settings.body_force_x)
    node.SetSolutionStepValue(BODY_FORCE_Y, 0, problem_settings.body_force_y)
    node.SetSolutionStepValue(BODY_FORCE_Z, 0, problem_settings.body_force_z)
    node.Free(PRESSURE)
    node.SetSolutionStepValue(PRESSURE, 0, 0.0)
    dist = node.GetSolutionStepValue(DISTANCE)
    if(dist < 0.0):
        n_active = n_active + 1
        node.SetSolutionStepValue(DISTANCE, 0, -small_value)
    else:
        node.SetSolutionStepValue(DISTANCE, 0, small_value)

if(n_active == 0):
    raise "ERROR. At least one node has to be initialized with a distance lesser than 0"

# make sure that the porosity is not zero on any node (set by default to fluid only)
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
# print body_force
viscosity = problem_settings.viscosity
density = problem_settings.density
fluid_solver = level_set_elembased_fluid_solver.ElemBasedLevelSetSolver(fluid_model_part, domain_size, body_force)

fluid_solver.redistance_frequency = problem_settings.redistance_frequency
fluid_solver.number_of_extrapolation_layers = problem_settings.extrapolation_layers

fluid_solver.Initialize()
#


print("fluid solver created")

# settings to be changed
# Dt = problem_settings.time_step
final_time = problem_settings.max_time
output_dt = problem_settings.output_dt
coef = problem_settings.delta_time_coefficient

# number_of_inital_steps = problem_settings.number_of_inital_steps
# initial_time_step = problem_settings.initial_time_step
out = 0


# mesh to be printed
if(problem_settings.print_layers == False):
    mesh_name = 0.0
    gid_io.InitializeMesh(mesh_name)
    gid_io.WriteMesh(fluid_model_part.GetMesh())
    gid_io.FinalizeMesh()
    gid_io.Flush()

    gid_io.InitializeResults(mesh_name, (fluid_model_part).GetMesh())

time_old_print = 0.0
time = 0.0
step = 0
initial_time_step = 0.00001
next_output_time = output_dt
Dt_old = problem_settings.time_step

while(time < final_time):

    if(step < 10):
        Dt = initial_time_step
    else:
        # Calculate Dt when a jump in velocity is reached
        Dt_new = fluid_solver.CalculateDelta_t(Dt)
        if(Dt_old >= coef * Dt_new):
            Dt = coef * Dt_new
        else:
            Dt = Dt_old

    time = time + Dt
    fluid_model_part.CloneTimeStep(time)

    print("******** CURRENT TIME = ", time)

    if(step >= 3):
# Calculate Dt when a jump in velocity is reached
# Dt_old = problem_settings.time_step
# Dt_new = fluid_solver.CalculateDelta_t(Dt)
# if(Dt_old >= coef * Dt_new):
# Dt = coef * Dt_new

        fluid_solver.Solve()

    time_to_print = time - time_old_print
#    if(time >= next_output_time):
    if(time_to_print >= problem_settings.output_dt):
        if(problem_settings.print_layers):
            # writing mesh
            gid_io.InitializeMesh(time)
            gid_io.WriteMesh((fluid_model_part).GetMesh())
            gid_io.FinalizeMesh()
            gid_io.InitializeResults(time, (fluid_model_part).GetMesh())

        gid_io.WriteNodalResults(PRESSURE, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(VELOCITY, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(DISTANCE, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(VISCOSITY, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(DENSITY, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(NORMAL, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(POROSITY, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(DIAMETER, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(IS_STRUCTURE, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(CONVECTION_VELOCITY, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(BODY_FORCE, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(AUX_INDEX, fluid_model_part.Nodes, time, 0)
        gid_io.Flush()

        if(problem_settings.print_layers):
            gid_io.FinalizeResults()
        time_old_print = time
#        next_output_time = time + output_dt

 #       out = 0
 #   out = out + 1
    step = step + 1

if(problem_settings.print_layers == False):
    gid_io.FinalizeResults()
