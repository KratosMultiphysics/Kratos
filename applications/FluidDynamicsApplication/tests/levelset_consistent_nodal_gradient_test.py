import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import math
import os

# from KratosMultiphysics.gid_output_process import GiDOutputProcess

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class ConsistentLevelsetNodalGradientTest(KratosUnittest.TestCase):

    def testConsistentGradientSquare2D(self):
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)

        KratosMultiphysics.ModelPartIO(GetFilePath("DistanceSmoothingTest/two_dim_symmetrical_square")).ReadModelPart(model_part)

        model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
        model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)
        model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 1.0)

        model_part.SetBufferSize(2)

        for node in model_part.Nodes:
            distance = math.sqrt((node.X+0.001)**2+(node.Y-0.001)**2) - 0.006
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, distance )
            node.SetValue(KratosMultiphysics.NODAL_AREA, 0.0)
            if (distance > 0.0):
                node.SetSolutionStepValue(KratosMultiphysics.PRESSURE,
                    2.0*node.X - 4.0*node.Y )
            else:
                node.SetSolutionStepValue(KratosMultiphysics.PRESSURE,
                    -4.0*node.X + 2.0*node.Y )

        KratosCFD.CalulateLevelsetConsistentNodalGradientProcess(model_part).Execute()

        node = (model_part.Nodes)[3]
        self.assertAlmostEqual(node.GetValue(KratosMultiphysics.PRESSURE_GRADIENT_X), -4.0)
        node = (model_part.Nodes)[7]
        self.assertAlmostEqual(node.GetValue(KratosMultiphysics.PRESSURE_GRADIENT_Y), 2.0)
        node = (model_part.Nodes)[30]
        self.assertAlmostEqual(node.GetValue(KratosMultiphysics.PRESSURE_GRADIENT_X), 2.0)
        node = (model_part.Nodes)[38]
        self.assertAlmostEqual(node.GetValue(KratosMultiphysics.PRESSURE_GRADIENT_Y), -4.0)

        # gid_output = GiDOutputProcess(model_part,
        #                            "consistent_gradient_test_2D",
        #                            KratosMultiphysics.Parameters("""
        #                                {
        #                                    "result_file_configuration" : {
        #                                        "gidpost_flags": {
        #                                            "GiDPostMode": "GiD_PostBinary",
        #                                            "WriteDeformedMeshFlag": "WriteUndeformed",
        #                                            "WriteConditionsFlag": "WriteConditions",
        #                                            "MultiFileFlag": "SingleFile"
        #                                        },
        #                                        "nodal_results"       : ["DISTANCE","PRESSURE"],
        #                                        "nodal_nonhistorical_results": ["PRESSURE_GRADIENT"]
        #                                    }
        #                                }
        #                                """)
        #                            )

        # gid_output.ExecuteInitialize()
        # gid_output.ExecuteBeforeSolutionLoop()
        # gid_output.ExecuteInitializeSolutionStep()
        # gid_output.PrintOutput()
        # gid_output.ExecuteFinalizeSolutionStep()
        # gid_output.ExecuteFinalize()

if __name__ == '__main__':
    KratosUnittest.main()