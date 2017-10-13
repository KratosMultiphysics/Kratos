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

solver.Dt = 1.0 / 32.0
solver.NumCells = [48, 48, 48]
solver.BorderWidth = [1, 1, 1]

solver.Initialize()

# Main loop
for i in range(0, 100):

    # Solve
    solver.SolveSolutionStep()

    # # Print results
    # if not i % 10:
    #     print("Printing results for step {}...".format(i))
    #     solver.WriteResults(i+1)
