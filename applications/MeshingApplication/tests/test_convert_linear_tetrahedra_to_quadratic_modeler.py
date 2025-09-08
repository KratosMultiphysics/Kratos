from importlib import import_module
from pathlib import Path

import KratosMultiphysics as kratos

from KratosMultiphysics import KratosUnittest
from KratosMultiphysics.MeshingApplication.modelers.convert_linear_tetrahedra_to_quadratic_modeler import ConvertLinearTetrahedraToQuadraticModeler


class TestConvertLinearTetrahedraToQuadraticModeler(KratosUnittest.TestCase):

    def setUp(self):
        self.model = kratos.Model()
        model_part = self.model.CreateModelPart("main_model_part")
        model_part.ProcessInfo.SetValue(kratos.DOMAIN_SIZE, 3)

        work_dir = Path(__file__).parent.absolute()
        mdpa_io = kratos.ModelPartIO(work_dir / "cube_with_5_faces")
        mdpa_io.ReadModelPart(model_part)


    def make_modeler(self, parameters: kratos.Parameters):
        registry_data = kratos.Registry[
            "Modelers.KratosMultiphysics.MeshingApplication.ConvertLinearTetrahedraToQuadraticModeler"
        ]
        module = import_module(registry_data["ModuleName"])
        return getattr(module, registry_data["ClassName"])(self.model, parameters)


    def test_creation_from_registry(self):
        modeler = self.make_modeler(kratos.Parameters('''{
            "model_part_name": "main_model_part"
        }'''))

        self.assertIsInstance(modeler, ConvertLinearTetrahedraToQuadraticModeler)


    def test_quadratic_refinement(self):

        model_part = self.model.GetModelPart("main_model_part")

        self.assertEqual(model_part.NumberOfNodes(), 64)
        self.assertEqual(model_part.NumberOfElements(), 162)
        self.assertEqual(model_part.NumberOfConditions(), 90)
        self.assertEqual(len(model_part.Elements[1].GetNodes()), 4)
        self.assertEqual(len(model_part.Conditions[163].GetNodes()), 3)

        modeler = ConvertLinearTetrahedraToQuadraticModeler(
            self.model,
            kratos.Parameters('''{
                "model_part_name": "main_model_part"
            }'''))

        modeler.SetupModelPart()

        self.assertEqual(model_part.NumberOfNodes(), 343)
        self.assertEqual(model_part.NumberOfElements(), 162)
        self.assertEqual(model_part.NumberOfConditions(), 90)
        self.assertEqual(len(model_part.Elements[1].GetNodes()), 10)
        self.assertEqual(len(model_part.Conditions[163].GetNodes()), 6)


    def test_root_required(self):

        modeler = ConvertLinearTetrahedraToQuadraticModeler(
            self.model,
            kratos.Parameters('''{
                "model_part_name": "main_model_part.MainPart.VolumeParts.Body1"
            }'''))

        msg = "A root ModelPart is required, got SubModelPart main_model_part.MainPart.VolumeParts.Body1"
        with self.assertRaisesRegex(Exception, msg):
            modeler.SetupModelPart()



if __name__ == "__main__":
    KratosUnittest.main()