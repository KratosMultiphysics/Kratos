from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
import math
import os

from KratosMultiphysics.gid_output_process import GiDOutputProcess

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestNodalDivergence(KratosUnittest.TestCase):
            
    def test_nodal_divergence_2d_square(self):
        print( " " )
        print( "Test 2D, square domain 1x1:" )
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        
        KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/two_dim_symmetrical_square")).ReadModelPart(model_part)
        
        model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
        model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)
        model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 1.0)

        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, (math.sqrt((node.X-0.5)**2+(node.Y-0.5)**2) - 0.5) )
            node.SetValue(KratosMultiphysics.NODAL_AREA, 0.0)

        KratosMultiphysics.ComputeNodalGradientProcess(
            model_part, 
            KratosMultiphysics.DISTANCE,
            KratosMultiphysics.DISTANCE_GRADIENT,
            KratosMultiphysics.NODAL_AREA).Execute()

        KratosMultiphysics.ComputeNodalNormalDivergenceProcess(
            model_part, 
            KratosMultiphysics.DISTANCE_GRADIENT,
            KratosMultiphysics.PRESSURE, #KratosCFD.CURVATURE,
            KratosMultiphysics.NODAL_AREA).Execute()

        for node in model_part.Nodes:
            if ( ((node.X-0.5)**2+(node.Y-0.5)**2 - 0.5**2) > -0.02 and ((node.X-0.5)**2+(node.Y-0.5)**2 - 0.5**2) < 0.02 ):
                print( node.GetSolutionStepValue(KratosMultiphysics.PRESSURE) )
            
        gid_output = GiDOutputProcess(model_part,
                                    "curvature_post_2D",
                                    KratosMultiphysics.Parameters("""
                                        {
                                            "result_file_configuration" : {
                                                "gidpost_flags": {
                                                    "GiDPostMode": "GiD_PostBinary",
                                                    "WriteDeformedMeshFlag": "WriteUndeformed",
                                                    "WriteConditionsFlag": "WriteConditions",
                                                    "MultiFileFlag": "SingleFile"
                                                },
                                                "nodal_results"       : ["PRESSURE"]
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
        
    def test_nodal_divergence_3d_cube(self):
        print( " " )
        print( "Test 3, cubic domain 1x1x1:" )
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        
        KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/three_dim_symmetrical_cube")).ReadModelPart(model_part)
        
        model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)
        model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 1.0)

        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, (math.sqrt((node.X-0.5)**2+(node.Y-0.5)**2+(node.Z-0.5)**2) - 0.5) )
            node.SetValue(KratosMultiphysics.NODAL_AREA, 0.0)

        KratosMultiphysics.ComputeNodalGradientProcess(
            model_part, 
            KratosMultiphysics.DISTANCE,
            KratosMultiphysics.DISTANCE_GRADIENT,
            KratosMultiphysics.NODAL_AREA).Execute()

        KratosMultiphysics.ComputeNodalNormalDivergenceProcess(
            model_part, 
            KratosMultiphysics.DISTANCE_GRADIENT,
            KratosMultiphysics.PRESSURE, #KratosCFD.CURVATURE,
            KratosMultiphysics.NODAL_AREA).Execute()
            
        for node in model_part.Nodes:
            if ( ((node.X-0.5)**2+(node.Y-0.5)**2+(node.Z-0.5)**2 - 0.5**2) > -0.01 and ((node.X-0.5)**2+(node.Y-0.5)**2+(node.Z-0.5)**2 - 0.5**2) < 0.01 ):
                print( node.GetSolutionStepValue(KratosMultiphysics.PRESSURE) )
                #self.assertAlmostEqual( node.GetSolutionStepValue(KratosMultiphysics.PRESSURE), 2.0/0.5, 1 )
            
        gid_output = GiDOutputProcess(model_part,
                                    "curvature_post_3D",
                                    KratosMultiphysics.Parameters("""
                                        {
                                            "result_file_configuration" : {
                                                "gidpost_flags": {
                                                    "GiDPostMode": "GiD_PostBinary",
                                                    "WriteDeformedMeshFlag": "WriteUndeformed",
                                                    "WriteConditionsFlag": "WriteConditions",
                                                    "MultiFileFlag": "SingleFile"
                                                },
                                                "nodal_results"       : ["PRESSURE"]
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

if __name__ == '__main__':
    KratosUnittest.main()
