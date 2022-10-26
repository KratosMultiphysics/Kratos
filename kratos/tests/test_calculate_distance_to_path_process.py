import KratosMultiphysics
import KratosMultiphysics.vtk_output_process as vtk_output_process

import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

current_model = KratosMultiphysics.Model()

# Analytic distance function
def DistanceFunction(node, radius):
    x = node.X
    y = node.Y
    node_radius = (x**2 + y**2)**0.5 - 1.0
    z = node.Z
    sub_radius = (node_radius**2 + z**2)**0.5 - radius
    return sub_radius - radius

# Import torus
model_part_torus = current_model.CreateModelPart("Torus")
model_part_torus.ProcessInfo[KratosMultiphysics.STEP] = 0
model_part_torus.CloneTimeStep(0.0)
model_part_torus.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
model_part_torus.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
input_mdpa = GetFilePath("torus3d")
model_part_io_torus = KratosMultiphysics.ModelPartIO(input_mdpa)
model_part_io_torus.ReadModelPart(model_part_torus)

# Import circle
radius = 0.1
model_part_circle = current_model.CreateModelPart("Circle")
model_part_circle.ProcessInfo[KratosMultiphysics.STEP] = 0
model_part_circle.CloneTimeStep(0.0)
model_part_circle.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
input_mdpa = GetFilePath("circle1d")
model_part_io_circle = KratosMultiphysics.ModelPartIO(input_mdpa)
model_part_io_circle.ReadModelPart(model_part_circle)

# Set the distance function
## Analytic distance function
for node in model_part_torus.Nodes:
    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, DistanceFunction(node, radius))

## Compute distance

# Output
vtk_output_parameters = KratosMultiphysics.Parameters("""{
    "Parameters" : {
        "model_part_name"                    : "Torus",
        "file_format"                        : "ascii",
        "output_precision"                   : 8,
        "output_interval"                    : 2,
        "output_sub_model_parts"             : false,
        "output_path"                        : "vtk_output_torus",
        "nodal_solution_step_data_variables" : ["DISTANCE", "TEMPERATURE"]
    }
}""")

vtk_output_process_torus = vtk_output_process.Factory(vtk_output_parameters, current_model)
vtk_output_process_torus.ExecuteInitializeSolutionStep()
vtk_output_process_torus.ExecuteFinalizeSolutionStep()
vtk_output_process_torus.PrintOutput()

vtk_output_parameters = KratosMultiphysics.Parameters("""{
    "Parameters" : {
        "model_part_name"                    : "Circle",
        "file_format"                        : "ascii",
        "output_precision"                   : 8,
        "output_interval"                    : 2,
        "output_sub_model_parts"             : false,
        "output_path"                        : "vtk_output_circle",
        "nodal_solution_step_data_variables" : ["TEMPERATURE"]
    }
}""")

vtk_output_process_circle = vtk_output_process.Factory(vtk_output_parameters, current_model)
vtk_output_process_circle.ExecuteInitializeSolutionStep()
vtk_output_process_circle.ExecuteFinalizeSolutionStep()
vtk_output_process_circle.PrintOutput()