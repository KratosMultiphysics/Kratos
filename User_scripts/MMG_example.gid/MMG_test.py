from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import os
import math
import KratosMultiphysics
import KratosMultiphysics.MeshingApplication as KratosMesh
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
from KratosMultiphysics.gid_output_process import GiDOutputProcess

class MMGParameterStudy():

    def __init__(self, Model_name, settings = KratosMultiphysics.Parameters("""{}""")):

        file_path = os.path.dirname(os.path.realpath(__file__))
        current_model = KratosMultiphysics.Model()
        self.model_part = current_model.CreateModelPart("MainModelPart")
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
        KratosMultiphysics.ModelPartIO(file_path+Model_name).ReadModelPart(self.model_part)
        #KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, self.model_part.Nodes)
        self.find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(self.model_part)
        self.local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(self.model_part, KratosMultiphysics.DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT, KratosMultiphysics.NODAL_AREA)
    
    def CalculateDistanceAndGradient(self):

        circle_radius = 20.0
        circle_center = [0.0, 0.0]
        for node in self.model_part.Nodes:
            distance = math.sqrt((node.X - circle_center[0])**2 + (node.Y - circle_center[1])**2) - circle_radius
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, distance)
        
        self.find_nodal_h.Execute()       
        self.local_gradient.Execute()
    
    def RemeshMetricProcess(self, h_min_max, min_size, bl_distance):

        ZeroVector = KratosMultiphysics.Vector(3)
        ZeroVector[0] = 0.0
        ZeroVector[1] = 0.0
        ZeroVector[2] = 0.0
        for node in self.model_part.Nodes:
            node.SetValue(KratosMesh.METRIC_TENSOR_2D, ZeroVector)

        # We define a metric using the ComputeLevelSetSolMetricProcess
        level_set_param = KratosMultiphysics.Parameters("""
                                {
                                    "minimal_size"                         : """+ str(min_size) +""",
                                    "enforce_current"                      : false,
                                    "anisotropy_remeshing"                 : true,
                                    "anisotropy_parameters":
                                    {
                                        "hmin_over_hmax_anisotropic_ratio"      : """+ str(h_min_max) +""",
                                        "boundary_layer_max_distance"           : """+ str(bl_distance) +""",
                                        "interpolation"                         : "Linear"
                                    }
                                }
                                """)
        metric_process = KratosMesh.ComputeLevelSetSolMetricProcess2D(self.model_part, KratosMultiphysics.DISTANCE_GRADIENT, level_set_param)
        metric_process.Execute()

        # We create the remeshing process
        remesh_param = KratosMultiphysics.Parameters("""{ }""")
        mmg_process = KratosMesh.MmgProcess2D(self.model_part, remesh_param)
        mmg_process.Execute()
    
    def CreateGIDOutput(self, model_name):

        gid_output = GiDOutputProcess(self.model_part,
                                    model_name,
                                    KratosMultiphysics.Parameters("""
                                        {
                                            "result_file_configuration" : {
                                                "gidpost_flags": {
                                                    "GiDPostMode": "GiD_PostBinary",
                                                    "WriteDeformedMeshFlag": "WriteUndeformed",
                                                    "WriteConditionsFlag": "WriteConditions",
                                                    "MultiFileFlag": "SingleFile"
                                                },
                                                "nodal_results"       : ["DISTANCE", "DISTANCE_GRADIENT"],
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
    
if __name__ == "__main__":

    h_min_max = [0.1, 0.05, 1.0]
    min_size = [0.25, 0.5, 0.75, 1.5, 2.5]
    bl_distance = [0.5, 0.25, 2.0]
    mmg_test_tool = MMGParameterStudy("/MMG_example")
    mmg_test_tool.CalculateDistanceAndGradient()
    mmg_test_tool.RemeshMetricProcess(h_min_max[2], min_size[0], bl_distance[2])
    mmg_test_tool.CalculateDistanceAndGradient()
    mmg_test_tool.CreateGIDOutput("output_05")

