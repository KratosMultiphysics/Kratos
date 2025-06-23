
import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.OptimizationApplication.model_part_controllers.mdpa_model_part_controller import MdpaModelPartController

class TestMdpaModelPartController(kratos_unittest.TestCase):
    def test_MdpaModelPartControllerReadData(self):
        with kratos_unittest.WorkFolderScope("model_part_utils_test", __file__):
            model_part_controller_settings = Kratos.Parameters("""{
                "model_part_name": "Structure",
                "input_filename" : "quads_with_nodal_values",
                "domain_size"    : 3,
                "read_data"      : true
            }""")

            model = Kratos.Model()
            model_part_controller = MdpaModelPartController(model, model_part_controller_settings)
            model_part_controller.GetModelPart().AddNodalSolutionStepVariable(Kratos.DENSITY)
            model_part_controller.ImportModelPart()

            ref_model_part = model.CreateModelPart("ref")
            Kratos.ModelPartIO("quads_with_nodal_values", Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(ref_model_part)

            self.assertEqual(model_part_controller.GetModelPart().NumberOfNodes(), ref_model_part.NumberOfNodes())
            self.assertEqual(model_part_controller.GetModelPart().NumberOfConditions(), ref_model_part.NumberOfConditions())
            self.assertEqual(model_part_controller.GetModelPart().NumberOfElements(), ref_model_part.NumberOfElements())

            node_id_value_map = {1: 10, 2: 20, 3: 30}
            for node in model_part_controller.GetModelPart().Nodes:
                if node.Id in node_id_value_map.keys():
                    self.assertEqual(node.GetSolutionStepValue(Kratos.DENSITY), node_id_value_map[node.Id])
                else:
                    self.assertEqual(node.GetSolutionStepValue(Kratos.DENSITY), 0.0)

    def test_MdpaModelPartController(self):
        with kratos_unittest.WorkFolderScope("model_part_utils_test", __file__):
            model_part_controller_settings = Kratos.Parameters("""{
                "model_part_name": "Structure",
                "input_filename" : "quads",
                "domain_size"    : 3
            }""")

            model = Kratos.Model()
            model_part_controller = MdpaModelPartController(model, model_part_controller_settings)
            model_part_controller.ImportModelPart()

            ref_model_part = model.CreateModelPart("ref")
            Kratos.ModelPartIO("quads", Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(ref_model_part)

            self.assertEqual(model_part_controller.GetModelPart().NumberOfNodes(), ref_model_part.NumberOfNodes())
            self.assertEqual(model_part_controller.GetModelPart().NumberOfConditions(), ref_model_part.NumberOfConditions())
            self.assertEqual(model_part_controller.GetModelPart().NumberOfElements(), ref_model_part.NumberOfElements())

    @classmethod
    def tearDownClass(cls) -> None:
        with kratos_unittest.WorkFolderScope("model_part_utils_test", __file__):
            DeleteFileIfExisting("quads.time")
            DeleteFileIfExisting("quads_with_nodal_values.time")

if __name__ == "__main__":
    kratos_unittest.main()