# Backward compatibility with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Python imports

# Kratos
from KratosMultiphysics import *
from KratosMultiphysics.FDApplication import *

# =========================================================================== #
# Procedures defined in this region should be moved to another file           #
# =========================================================================== #


def GenerateGidIO(params):

    # Create a dictionary with all the values for the gid_io config
    gid_mode = {
        "Binary": GiDPostMode.GiD_PostBinary,
        "Ascii":  GiDPostMode.GiD_PostAscii,
        "AsciiZipped": GiDPostMode.GiD_PostAsciiZipped,
    }

    multifile = {
        "Single": MultiFileFlag.SingleFile,
        "Multiples": MultiFileFlag.MultipleFiles,
    }

    deformed_mesh_flag = {
        "Deformed": WriteDeformedMeshFlag.WriteDeformed,
        "Undeformed": WriteDeformedMeshFlag.WriteUndeformed,
    }

    write_conditions = {
        "OnlyElements": WriteConditionsFlag.WriteElementsOnly,
        "Conditions": WriteConditionsFlag.WriteConditions,
    }

    gid_io = GidIO(
        params["FileName"],
        gid_mode[params["GiDMode"]],
        multifile[params["MultiFileFlag"]],
        deformed_mesh_flag[params["WriteDeformedMeshFlag"]],
        write_conditions[params["WriteConditionsFlag"]],
    )

    return gid_io


# =========================================================================== #
#                                                                             #
# =========================================================================== #

# Define the varianles to be used
variables_dictionary = {
    "PRESSURE": PRESSURE,
    "VELOCITY": VELOCITY,
}

# Modelpart for the fluid
model_part = ModelPart("FluidPart")

model_part.AddNodalSolutionStepVariable(PRESSURE)
model_part.AddNodalSolutionStepVariable(VELOCITY)
model_part.AddNodalSolutionStepVariable(VISCOSITY)
model_part.AddNodalSolutionStepVariable(DENSITY)

# initialize GiD  I/O
input_file_name = "IsoCubeFluid"  # ProjectParameters.problem_name

# reading the fluid part
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(model_part)

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
model_part.SetBufferSize(1)

# Configure GiD  I/O
gid_io_configuration = {
    "FileName": input_file_name,
    "GiDMode": "Binary",              # Binary, Ascii, AsciiZipped
    "WriteDeformedMeshFlag": "algo",  # WriteDeformed. WriteUndeformed
    "WriteConditionsFlag": "algo2",   # WriteConditions, WriteElementsOnly
    "MultiFileFlag": "algo3",         # MultiFileFlag.SingleFile MultipleFiles.
}

# gid_io = GenerateGidIO(gid_io_configuration)

solver = BfeccSolverStrategy(model_part)

solver.Initialize()

for i in range(0, 500):
    # print("Convecting iteration {}".format(i))
    solver.SolveSolutionStep()

# # Stepping and time settings
# Dt = ProjectParameters.Dt
# full_Dt = Dt
# initial_Dt = 0.001 * full_Dt  # 0.05 #0.01
# Nsteps = ProjectParameters.nsteps
# final_time = ProjectParameters.max_time
# output_time = ProjectParameters.output_time
#
# time = ProjectParameters.Start_time
# out = 0
# step = 0
#
# while(time <= final_time):
#
#     if(step < 3):
#         Dt = initial_Dt
#     else:
#         Dt = full_Dt
#
#     time = time + Dt
#     step = step + 1
#     fluid_model_part.CloneTimeStep(time)
#
#     print("STEP = ", step)
#     print("TIME = ", time)
#
#     if(step >= 3):
#         fluid_solver.Solve()
#
#         if(step < 4):
#             for k in range(0, ProjectParameters.divergence_cleareance_step):
#                 print("DOING DIVERGENCE CLEAREANCE")
#                 buffer_size = fluid_model_part.GetBufferSize()
#                 for i in range(0, buffer_size):
#                     for node in fluid_model_part.Nodes:
#                         vel = node.GetSolutionStepValue(VELOCITY)
#                         node.SetSolutionStepValue(VELOCITY, i, vel)
#                         node.SetSolutionStepValue(PRESSURE, i, 0.0)
#                     if(SolverType == "monolithic_solver_eulerian"):
#                         zero_vector = Vector(3)
#                         zero_vector[0] = 0.0
#                         zero_vector[1] = 0.0
#                         zero_vector[2] = 0.0
#                         for node in fluid_model_part.Nodes:
#                             node.SetSolutionStepValue(ACCELERATION, i, zero_vector)
#                     if(ProjectParameters.TurbulenceModel == "Spalart-Allmaras"):
#                         for node in fluid_model_part.Nodes:
#                             visc = node.GetSolutionStepValue(VISCOSITY)
#                             node.SetSolutionStepValue(VISCOSITY, i, visc)
#
#                 fluid_solver.Solve()
#
#         graph_printer.PrintGraphs(time)
#         PrintDrag(drag_list, drag_file_output_list, fluid_model_part, time)
#
#     if(output_time <= out):
#         if(ProjectParameters.VolumeOutput):
#             if Multifile:
#                 gid_io.InitializeMesh(time)
#                 gid_io.WriteMesh(fluid_model_part.GetMesh())
#                 gid_io.FinalizeMesh()
#
#                 gid_io.InitializeResults(time, (fluid_model_part).GetMesh())
#
#             PrintResults(fluid_model_part)
#             out = 0
#         else:
#             cut_model_part.CloneTimeStep(time)
#             Cut_App.UpdateCutData(cut_model_part, fluid_model_part)
#
#             if Multifile:
#                 gid_io.InitializeMesh(time)
#                 gid_io.WriteMesh(cut_model_part.GetMesh())
#                 gid_io.FinalizeMesh()
#
#                 gid_io.InitializeResults(time, (cut_model_part).GetMesh())
#
#             PrintResults(cut_model_part)
#
#             out = 0
#
#         if Multifile:
#             gid_io.FinalizeResults()
#             if ProjectParameters.GiDPostMode == "Binary":
#                 f.write(ProjectParameters.problem_name + '_' + str(time) + '.post.bin\n')
#             elif ProjectParameters.GiDPostMode == "Ascii":
#                 f.write(ProjectParameters.problem_name + '_' + str(time) + '.post.msh\n')
#
#     out = out + Dt
#
# if Multifile:
#     f.close()
# else:
#     gid_io.FinalizeResults()
#
# for i in drag_file_output_list:
#     i.close();
