from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#
#
# import the configuration data as read from the GiD
import ProjectParameters
import define_output


#
#
# setting the domain size for the problem to be solved
domain_size = ProjectParameters.domain_size

#
#
import sys
sys.path.append(ProjectParameters.kratos_path)
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *

# defining variables to be used
# GID IO IS NOT USING THIS NOW. TO BE REMOVED ONCE THE "PRINT IN POINTS"
# CODE IS NOT USING IT

variables_dictionary = {"PRESSURE": PRESSURE,
                        "VELOCITY": VELOCITY,
                        "REACTION": REACTION,
                        "DISTANCE": DISTANCE, }

# defining a model part for the fluid
fluid_model_part = ModelPart("FluidPart")

if "REACTION" in ProjectParameters.nodal_results:
    fluid_model_part.AddNodalSolutionStepVariable(REACTION)
if "DISTANCE" in ProjectParameters.nodal_results:
    fluid_model_part.AddNodalSolutionStepVariable(DISTANCE)

#
#
# importing the solvers needed
SolverSettings = ProjectParameters.FluidSolverConfiguration
solver_module = import_solver(SolverSettings)

#
#
# importing variables
solver_module.AddVariables(fluid_model_part, SolverSettings)

# introducing input file name
input_file_name = ProjectParameters.problem_name

# reading the fluid part
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(3)

# Check tetrahedral mesh for wrong orientation
throw_errors = False
orientation_check = TetrahedralMeshOrientationCheck(fluid_model_part,throw_errors)
orientation_check.Execute()

solver_module.AddDofs(fluid_model_part, SolverSettings)

# If Lalplacian form = 2, free all pressure Dofs
# laplacian_form = ProjectParameters.laplacian_form
# if(laplacian_form >= 2):
    # for node in fluid_model_part.Nodes:
        # node.Free(PRESSURE)

# copy Y_WALL
for node in fluid_model_part.Nodes:
    y = node.GetSolutionStepValue(Y_WALL, 0)
    node.SetValue(Y_WALL, y)

#
#
# Creating the fluid solver
fluid_solver = solver_module.CreateSolver(
    fluid_model_part, SolverSettings)

# activate turbulence model
if(ProjectParameters.FluidSolverConfiguration.TurbulenceModel == "Spalart-Allmaras"):
    # apply the initial turbulent viscosity on all of the nodes
    turb_visc = SolverSettings.TurbulentViscosity
    for node in fluid_model_part.Nodes:
        node.SetSolutionStepValue(TURBULENT_VISCOSITY, 0, turb_visc)
        visc = node.GetSolutionStepValue(VISCOSITY)
        node.SetSolutionStepValue(MOLECULAR_VISCOSITY, 0, visc)
        if (node.IsFixed(VELOCITY_X) and node.GetSolutionStepValue(VELOCITY_X, 0) != 0.0) or \
           (node.IsFixed(VELOCITY_Y) and node.GetSolutionStepValue(VELOCITY_Y, 0) != 0.0) or \
           (node.IsFixed(VELOCITY_Z) and node.GetSolutionStepValue(VELOCITY_Z, 0) != 0.0):
            node.Fix(TURBULENT_VISCOSITY)

    # select nodes on the wall
    fluid_solver.wall_nodes = []
    for i in SolverSettings.SA_wall_group_ids:
        # get the nodes of the wall for SA.
        nodes = fluid_model_part.GetNodes(i)
        for node in nodes:
            fluid_solver.wall_nodes.append(node)
            node.SetSolutionStepValue(TURBULENT_VISCOSITY, 0, 0.0)
            node.Fix(TURBULENT_VISCOSITY)


fluid_solver.Initialize()
print("fluid solver created")
#
#

# initialize GiD  I/O
from gid_output import GiDOutput
gid_io = GiDOutput(input_file_name,
                   ProjectParameters.VolumeOutput,
                   ProjectParameters.GiDPostMode,
                   ProjectParameters.GiDMultiFileFlag,
                   ProjectParameters.GiDWriteMeshFlag,
                   ProjectParameters.GiDWriteConditionsFlag)

if not ProjectParameters.VolumeOutput:
    cut_list = define_output.DefineCutPlanes()
    gid_io.define_cuts(fluid_model_part, cut_list)

gid_io.initialize_results(fluid_model_part)

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
    tmp = "#time RX RY RZ\n"
    f.write(tmp)
    f.flush()

print(drag_file_output_list)


def PrintDrag(drag_list, drag_file_output_list, fluid_model_part, time):
    i = 0
    for it in drag_list:
        print(it[0])
        nodes = fluid_model_part.GetNodes(it[0])
        drag = Vector(3)
        drag[0] = 0.0
        drag[1] = 0.0
        drag[2] = 0.0
        for node in nodes:
            reaction = node.GetSolutionStepValue(REACTION, 0)
            drag[0] += reaction[0]
            drag[1] += reaction[1]
            drag[2] += reaction[2]

        output = str(time) + " " + str(drag[0]) + " " + str(
            drag[1]) + " " + str(drag[2]) + "\n"
        # print drag_file_output_list[i]
        # print output
        drag_file_output_list[i].write(output)
        drag_file_output_list[i].flush()
        i = i + 1


# 33
# preparing output of point graphs
import point_graph_printer

output_nodes_list = define_output.DefineOutputPoints()
graph_printer = point_graph_printer.PrintGraphPrinter(
    output_nodes_list,
    fluid_model_part,
    variables_dictionary,
    domain_size)


# Stepping and time settings
Dt = ProjectParameters.Dt
Nsteps = ProjectParameters.nsteps
final_time = ProjectParameters.max_time
output_time = ProjectParameters.output_time

time = ProjectParameters.Start_time
out = 0
step = 0

while(time <= final_time):

    time = time + Dt
    step = step + 1
    fluid_model_part.CloneTimeStep(time)

    print("STEP = ", step)
    print("TIME = ", time)

    if(step >= 3):
        fluid_solver.Solve()

        graph_printer.PrintGraphs(time)
        PrintDrag(drag_list, drag_file_output_list, fluid_model_part, time)

    if(output_time <= out):
        gid_io.write_results(
            time,
            fluid_model_part,
            ProjectParameters.nodal_results,
            ProjectParameters.gauss_points_results)
        out = 0

    out = out + Dt

gid_io.finalize_results()

for i in drag_file_output_list:
    i.close()
