# We import the libraies
from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.MeshingApplication as MeshingApplication

import KratosMultiphysics.KratosUnittest as KratosUnittest

import os

class TestRemeshMMG(KratosUnittest.TestCase):
    
    def test_remesh_sphere(self):
        # We create the model part
        main_model_part = KratosMultiphysics.ModelPart("MainModelPart")
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)

        # We add the variables needed 
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)

        for node in main_model_part.Nodes:
            node.AddDof(KratosMultiphysics.DISTANCE)

        # We import the model main_model_part
        file_path = os.path.dirname(os.path.realpath(__file__))
        KratosMultiphysics.ModelPartIO(file_path + "/mmg_eulerian_test/coarse_sphere_test").ReadModelPart(main_model_part)

        # We calculate the gradient of the distance variable
        find_nodal_h = KratosMultiphysics.FindNodalHProcess(main_model_part)
        find_nodal_h.Execute()
        local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(main_model_part, KratosMultiphysics.DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT, KratosMultiphysics.NODAL_AREA)
        local_gradient.Execute()

        # We set to zero the metric
        ZeroVector = KratosMultiphysics.Vector(6) 
        ZeroVector[0] = 0.0
        ZeroVector[1] = 0.0
        ZeroVector[2] = 0.0
        ZeroVector[3] = 0.0
        ZeroVector[4] = 0.0
        ZeroVector[5] = 0.0

        for node in main_model_part.Nodes:
            node.SetValue(MeshingApplication.MMG_METRIC, ZeroVector)
                
        # We define a metric using the ComputeLevelSetSolMetricProcess
        MetricParameters = KratosMultiphysics.Parameters("""
        {
            "minimal_size"                      : 1.0e-1,
            "enforce_current"                   : false,
            "anisotropy_remeshing"              : false,
            "anisotropy_parameters"             :{
                "hmin_over_hmax_anisotropic_ratio"  : 0.15,
                "boundary_layer_max_distance"       : 1.0e-4,
                "interpolation"                     : "Linear"
            }
        }
        """)
        metric_process = MeshingApplication.ComputeLevelSetSolMetricProcess3D(main_model_part, KratosMultiphysics.DISTANCE_GRADIENT, MetricParameters)

        metric_process.Execute()

        mmg_parameters = KratosMultiphysics.Parameters("""
        {
            "filename"                         : "mmg_eulerian_test/coarse_sphere_test", 
            "save_external_files"              : true,
            "echo_level"                       : 0
        }
        """)

        # We create the remeshing utility
        mmg_parameters["filename"].SetString(file_path + "/" + mmg_parameters["filename"].GetString())
        mmg_process = MeshingApplication.MmgProcess3D(main_model_part, mmg_parameters)

        # We remesh
        mmg_process.Execute()
        
        # Finally we export to GiD
        from gid_output_process import GiDOutputProcess
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
                                                "nodal_results"       : []
                                            }
                                        }
                                        """)
                                    )

        #gid_output.ExecuteInitialize()
        #gid_output.ExecuteBeforeSolutionLoop()
        #gid_output.ExecuteInitializeSolutionStep()
        #gid_output.PrintOutput()
        #gid_output.ExecuteFinalizeSolutionStep()
        #gid_output.ExecuteFinalize()  
        
        import filecmp 
        value = filecmp.cmp(file_path + "/mmg_eulerian_test/coarse_sphere_test_result.mesh", file_path + "/mmg_eulerian_test/coarse_sphere_test_step=0.o.mesh")
        self.assertTrue(value)
        
if __name__ == '__main__':
    KratosUnittest.main()
