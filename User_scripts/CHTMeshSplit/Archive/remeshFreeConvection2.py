from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.MeshingApplication as MeshingApplication
import KratosMultiphysics.FluidDynamicsApplication as CFDApplication
import math

current_model = KratosMultiphysics.Model()
main_model_part = current_model.CreateModelPart("MainModelPart")
main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)

# We add the variables needed
main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)

# We import the model main_model_part
KratosMultiphysics.ModelPartIO("cavity0_2").ReadModelPart(main_model_part)
#KratosMultiphysics.ModelPartIO("test_embedded_2D_coarse").ReadModelPart(main_model_part)

find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(main_model_part)
find_nodal_h.Execute()

## Circle distance function
circle_radious = 10
center_coordinates = [0.0, 0.0]

for node in main_model_part.Nodes:
    distance = ((node.X-center_coordinates[0])**2+(node.Y-center_coordinates[1])**2)**0.5 - circle_radious
    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, distance)

# We calculate the gradient of the distance variable
local_gradient = KratosMultiphysics.ComputeNodalGradientProcess2D(main_model_part, KratosMultiphysics.DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT, KratosMultiphysics.NODAL_AREA)
local_gradient.Execute()

# We set to zero the metric
ZeroVector = KratosMultiphysics.Vector(3)
ZeroVector[0] = 0.0
ZeroVector[1] = 0.0
ZeroVector[2] = 0.0
for node in main_model_part.Nodes:
    node.SetValue(MeshingApplication.METRIC_TENSOR_2D, ZeroVector)

# We define a metric using the ComputeLevelSetSolMetricProcess
level_set_param = KratosMultiphysics.Parameters("""
                        {
                            "minimal_size"                         : 0.8,
                            "enforce_current"                      : false,
                            "anisotropy_remeshing"                 : true,
                            "anisotropy_parameters":
                            {
                                "hmin_over_hmax_anisotropic_ratio"      : 0.1,
                                "boundary_layer_max_distance"           : 1,
                                "interpolation"                         : "Linear"
                            }
                        }
                        """)
metric_process = MeshingApplication.ComputeLevelSetSolMetricProcess2D(main_model_part,KratosMultiphysics.DISTANCE_GRADIENT,level_set_param)
metric_process.Execute()

# We create the remeshing process
remesh_param = KratosMultiphysics.Parameters("""{ }""")
mmg_process = MeshingApplication.MmgProcess2D(main_model_part, remesh_param)
mmg_process.Execute()

for node in main_model_part.Nodes:
    distance = ((node.X-center_coordinates[0])**2+(node.Y-center_coordinates[1])**2)**0.5 - circle_radious
    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, distance)

# We calculate the gradient of the distance variable
local_gradient.Execute()

# Finally we export to GiD
from KratosMultiphysics.gid_output_process import GiDOutputProcess
gid_output = GiDOutputProcess(main_model_part,
                            "gid_output",
                            KratosMultiphysics.Parameters("""
                                {
                                    "result_file_configuration" : {
                                        "gidpost_flags": {
                                            "GiDPostMode": "GiD_PostBinary",
                                            "WriteDeformedMeshFlag": "WriteUndeformed",
                                            "WriteConditionsFlag": "WriteConditions",
                                            "MultiFileFlag": "SingleFile"
                                        },
                                        "nodal_results"       : ["DISTANCE","DISTANCE_GRADIENT"],
                                        "nodal_nonhistorical_results"       : ["ANISOTROPIC_RATIO"]
                                    }
                                }
                                """)
                            )

gid_output.ExecuteInitialize()
gid_output.ExecuteBeforeSolutionLoop()
gid_output.ExecuteInitializeSolutionStep()
gid_output.PrintOutput()
gid_output.ExecuteFinalizeSolutionStep()
gid_output.ExecuteFinalize()

# Finally we export to VTK
vtk_settings = KratosMultiphysics.Parameters("""{
    "model_part_name"                    : "MainModelPart_1",
    "file_format"                        : "ascii",
    "output_precision"                   : 7,
    "output_control_type"                : "step",
    "output_frequency"                   : 1.0,
    "output_sub_model_parts"             : false,
    "save_output_files_in_folder"        : false,
    "nodal_solution_step_data_variables" : ["DISTANCE","DISTANCE_GRADIENT"],
    "nodal_data_value_variables"         : ["ANISOTROPIC_RATIO"],
    "element_data_value_variables"       : [],
    "condition_data_value_variables"     : [],
    "gauss_point_variables"              : []
}""")

vtk_io = KratosMultiphysics.VtkOutput(main_model_part, vtk_settings)
vtk_io.PrintOutput()
