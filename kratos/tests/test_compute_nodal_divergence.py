import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
import math
import os

#from KratosMultiphysics.gid_output_process import GiDOutputProcess

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestNodalDivergence(KratosUnittest.TestCase):

    def test_nodal_normal_divergence_2d_square(self):
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
            KratosMultiphysics.PRESSURE,
            KratosMultiphysics.NODAL_AREA).Execute()

        node = (model_part.Nodes)[19]
        self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.PRESSURE), 2.1431154720650256)
        node = (model_part.Nodes)[18]
        self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.PRESSURE), 2.143115472065025)
        node = (model_part.Nodes)[9]
        self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.PRESSURE), 2.1431154720650243)

        #gid_output = GiDOutputProcess(model_part,
        #                            "curvature_post_2D",
        #                            KratosMultiphysics.Parameters("""
        #                                {
        #                                    "result_file_configuration" : {
        #                                        "gidpost_flags": {
        #                                            "GiDPostMode": "GiD_PostBinary",
        #                                            "WriteDeformedMeshFlag": "WriteUndeformed",
        #                                            "WriteConditionsFlag": "WriteConditions",
        #                                           "MultiFileFlag": "SingleFile"
        #                                        },
        #                                        "nodal_results"       : ["PRESSURE"]
        #                                    }
        #                                }
        #                                """)
        #                            )

        #gid_output.ExecuteInitialize()
        #gid_output.ExecuteBeforeSolutionLoop()
        #gid_output.ExecuteInitializeSolutionStep()
        #gid_output.PrintOutput()
        #gid_output.ExecuteFinalizeSolutionStep()
        #gid_output.ExecuteFinalize()

    def test_nodal_normal_divergence_3d_cube(self):
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

        node = (model_part.Nodes)[117]
        self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.PRESSURE), 3.7011835866693046)
        node = (model_part.Nodes)[36]
        self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.PRESSURE), 3.9753094244054954)
        node = (model_part.Nodes)[175]
        self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.PRESSURE), 3.701183586669305)

        #gid_output = GiDOutputProcess(model_part,
        #                            "curvature_post_3D",
        #                            KratosMultiphysics.Parameters("""
        #                                {
        #                                    "result_file_configuration" : {
        #                                        "gidpost_flags": {
        #                                            "GiDPostMode": "GiD_PostBinary",
        #                                            "WriteDeformedMeshFlag": "WriteUndeformed",
        #                                            "WriteConditionsFlag": "WriteConditions",
        #                                            "MultiFileFlag": "SingleFile"
        #                                        },
        #                                        "nodal_results"       : ["PRESSURE"]
        #                                    }
        #                                }
        #                                """)
        #                            )

        #gid_output.ExecuteInitialize()
        #gid_output.ExecuteBeforeSolutionLoop()
        #gid_output.ExecuteInitializeSolutionStep()
        #gid_output.PrintOutput()
        #gid_output.ExecuteFinalizeSolutionStep()
        #gid_output.ExecuteFinalize()

    def test_nodal_divergence_2d_square(self):
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DETERMINANT_F)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)

        KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/two_dim_symmetrical_square")).ReadModelPart(model_part)

        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, [(node.X-0.5)**2,0,0] )
            node.SetValue(KratosMultiphysics.NODAL_AREA, 0.0)

        normalize_vector_field = False
        KratosMultiphysics.ComputeNodalNormalDivergenceProcess(
            model_part,
            KratosMultiphysics.DISPLACEMENT,
            KratosMultiphysics.DETERMINANT_F,
            KratosMultiphysics.NODAL_AREA,
            normalize_vector_field).Execute()

        node = (model_part.Nodes)[19]
        self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DETERMINANT_F), 0.2)
        node = (model_part.Nodes)[18]
        self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DETERMINANT_F), -0.85)
        node = (model_part.Nodes)[9]
        self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DETERMINANT_F), -0.2)

    def test_nodal_divergence_3d_cube(self):
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DETERMINANT_F)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)

        KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/three_dim_symmetrical_cube")).ReadModelPart(model_part)

        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, [0,0,(node.Z-0.5)**2] )
            node.SetValue(KratosMultiphysics.NODAL_AREA, 0.0)

        normalize_vector_field = False
        KratosMultiphysics.ComputeNodalNormalDivergenceProcess(
            model_part,
            KratosMultiphysics.DISPLACEMENT,
            KratosMultiphysics.DETERMINANT_F,
            KratosMultiphysics.NODAL_AREA,
            normalize_vector_field).Execute()

        node = (model_part.Nodes)[117]
        self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DETERMINANT_F), 0.666666666)
        node = (model_part.Nodes)[36]
        self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DETERMINANT_F), 0.0)
        node = (model_part.Nodes)[175]
        self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DETERMINANT_F),-0.666666666)

if __name__ == '__main__':
    KratosUnittest.main()
