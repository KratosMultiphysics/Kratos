#
#
# setting the domain size for the problem to be solved
domain_size = 2

#
#
# ATTENTION: here the order is important

from KratosMultiphysics import *
from KratosMultiphysics.ConvectionDiffusionApplication import *

# defining a model part
model_part = ModelPart("FluidPart")

# read solver settings from configuration file
import ProjectParameters
SolverSettings = ProjectParameters.SolverSettings2

# import solver file
solver_constructor = __import__(SolverSettings.solver_type)

# define variables to be stored
solver_constructor.AddVariables(model_part, SolverSettings)

# introducing input file name
input_file_name = "face_heat"
# reading the fluid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO(input_file_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(model_part)

mesh_name = 0.0
gid_io.InitializeMesh(mesh_name)
gid_io.WriteMesh((model_part).GetMesh())
gid_io.FinalizeMesh()
print(model_part)

# the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(3)

# define dofs to be stored
solver_constructor.AddDofs(model_part, SolverSettings)

# construct the solver
conv_diff_solver = solver_constructor.CreateSolver(model_part, SolverSettings)
conv_diff_solver.Initialize()


# assigning the fluid properties
conductivity = 25.0
density = 900.0
specific_heat = 2400.0
temperature = 298.0
for node in model_part.Nodes:
    node.SetSolutionStepValue(CONDUCTIVITY, 0, conductivity);
    node.SetSolutionStepValue(DENSITY, 0, density);
    node.SetSolutionStepValue(SPECIFIC_HEAT, 0, specific_heat);
    node.SetSolutionStepValue(TEMPERATURE, 0, temperature);
    node.SetSolutionStepValue(VELOCITY_X, 0, 0.0);
    node.SetSolutionStepValue(VELOCITY_Y, 0, 0.0);


model_part.Properties[0][EMISSIVITY] = 0.0
model_part.Properties[0][AMBIENT_TEMPERATURE] = 0.0
model_part.Properties[0][CONVECTION_COEFFICIENT] = 0.0

# applying a temperature of 100
for node in model_part.Nodes:
    if(node.Y > 0.499):
        node.SetSolutionStepValue(FACE_HEAT_FLUX, 0, 1000.0);

# settings to be changed
Dt = ProjectParameters.Dt
full_Dt = Dt
initial_Dt = 0.01 * full_Dt  # 0.05 #0.01
Nsteps = ProjectParameters.nsteps
final_time = ProjectParameters.max_time
output_time = ProjectParameters.output_time
output_step = ProjectParameters.output_step
time = ProjectParameters.Start_time

out = 0
step = 0
output_step = 1
Nsteps = 3000
for step in range(0, Nsteps):
    time = Dt * step
    model_part.CloneTimeStep(time)

    # solving the fluid problem
    if(step > 3):
        conv_diff_solver.Solve()

    # print the results
    if(out == output_step):
        gid_io.WriteNodalResults(TEMPERATURE, model_part.Nodes, time, 0)
        # gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(TEMP_CONV_PROJ, model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(FACE_HEAT_FLUX, model_part.Nodes, time, 0)
        out = 0
    out = out + 1
