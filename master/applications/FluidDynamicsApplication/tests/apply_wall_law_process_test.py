# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
from KratosMultiphysics.FluidDynamicsApplication import apply_wall_law_process

class ApplyWallLawProcessTest(KratosUnittest.TestCase):
    @classmethod
    def __CreateModel(self):
        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart("TestModelPart")
        properties = model_part.CreateNewProperties(0)
        model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 2
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        wall_model_part = model_part.CreateSubModelPart("WallModelPart")
        wall_model_part.CreateNewNode(1, 0.0, 1.0, 0.0)
        wall_model_part.CreateNewNode(2, 2.0, 0.0, 0.0)
        wall_model_part.CreateNewCondition("LineCondition2D2N", 1, [1, 2], properties)
        return model

    @classmethod
    def __GetBlankParameters(cls):
        return KratosMultiphysics.Parameters("""{
            "python_module" : "apply_wall_law_process",
            "kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
            "process_name" : "ApplyWallLawProcess",
            "Parameters"    : {
                "model_part_name" : "TestModelPart.WallModelPart",
                "wall_model_name" : "",
                "wall_model_settings" : {}
            }
        }""")

    def testApplyWallLawProcessNavierSlip(self):
        model = self.__CreateModel()
        settings = self.__GetBlankParameters()
        settings["Parameters"]["wall_model_name"].SetString("navier_slip")
        settings["Parameters"]["wall_model_settings"].AddEmptyValue("slip_length").SetDouble(1.0)

        process = apply_wall_law_process.Factory(settings, model)
        process.ExecuteInitialize()

        wall_model_part = model.GetModelPart("TestModelPart.WallModelPart")
        self.assertEqual(wall_model_part.NumberOfConditions(), 1)
        for node in wall_model_part.Nodes:
            self.assertAlmostEqual(node.GetValue(KratosCFD.SLIP_LENGTH), 1.0)

    def testApplyWallLawProcessLinearLog(self):
        model = self.__CreateModel()
        settings = self.__GetBlankParameters()
        settings["Parameters"]["wall_model_name"].SetString("linear_log")
        settings["Parameters"]["wall_model_settings"].AddEmptyValue("y_wall").SetDouble(1.0)

        process = apply_wall_law_process.Factory(settings, model)
        process.ExecuteInitialize()

        wall_model_part = model.GetModelPart("TestModelPart.WallModelPart")
        self.assertEqual(wall_model_part.NumberOfConditions(), 1)
        for node in wall_model_part.Nodes:
            self.assertAlmostEqual(node.GetValue(KratosCFD.Y_WALL), 1.0)

if __name__ == '__main__':
    KratosUnittest.main()
