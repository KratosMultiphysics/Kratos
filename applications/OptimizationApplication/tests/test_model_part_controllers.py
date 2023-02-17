
import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest
import KratosMultiphysics.StructuralMechanicsApplication
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.model_part_controllers.mdpa_model_part_controller import MdpaModelPartController

@kratos_unittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
class TestMdpaModelPartController(kratos_unittest.TestCase):
    def test_MdpaModelPartController(self):
        with kratos_unittest.WorkFolderScope("linear_element_test", __file__):
            model_part_controller_settings = Kratos.Parameters("""{
                "model_part_name": "Structure",
                "input_filename" : "Structure",
                "domain_size"    : 3
            }""")

            model = Kratos.Model()
            model_part_controller = MdpaModelPartController(model, model_part_controller_settings, OptimizationInfo())
            model_part_controller.ImportModelPart()

            ref_model_part = model.CreateModelPart("ref")
            Kratos.ModelPartIO("Structure", Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(ref_model_part)

            self.assertEqual(model_part_controller.GetModelPart().NumberOfNodes(), ref_model_part.NumberOfNodes())
            self.assertEqual(model_part_controller.GetModelPart().NumberOfConditions(), ref_model_part.NumberOfConditions())
            self.assertEqual(model_part_controller.GetModelPart().NumberOfElements(), ref_model_part.NumberOfElements())

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()