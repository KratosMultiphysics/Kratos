
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.StructuralMechanicsApplication as KratosSMA

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.controls.shape.vertex_morphing_shape_control import VertexMorphingShapeControl

class TestShellShapeControl(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("shell")
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

        parameters = Kratos.Parameters("""{
            "controlled_model_part_names": ["shell.design"],
            "filter_settings": {
                "type" : "surface_implicit",
                "radius": 0.2
            },
            "fixed_model_parts" : {
                "shell.top_edge" : [true,true,true],
                "shell.edge_support" : [true,true,true]
            }
        }""")

        cls.optimization_problem = OptimizationProblem()
        cls.implicit_shape_control = VertexMorphingShapeControl("implicit", cls.model, parameters, cls.optimization_problem)
        cls.optimization_problem.AddComponent(cls.implicit_shape_control)
        parameters["filter_settings"]["type"].SetString("surface_explicit")
        cls.explicit_shape_control = VertexMorphingShapeControl("explicit", cls.model, parameters, cls.optimization_problem)
        cls.optimization_problem.AddComponent(cls.explicit_shape_control)

        with kratos_unittest.WorkFolderScope(".", __file__):
            Kratos.ModelPartIO("shell", Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(cls.model_part)
            material_settings = Kratos.Parameters("""{"Parameters": {"materials_filename": "materials_2D.json"}} """)
            Kratos.ReadMaterialsUtility(material_settings, cls.model)

        cls.implicit_shape_control.Initialize()
        cls.explicit_shape_control.Initialize()

    @classmethod
    def tearDownClass(cls):
        with kratos_unittest.WorkFolderScope(".", __file__):
            DeleteFileIfExisting("shell.time")

    def test_GetControlField(self):
        implicit_control_field = self.implicit_shape_control.GetControlField()
        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormInf(implicit_control_field), 1.0903407171261381, 10)
        explicit_control_field = self.explicit_shape_control.GetControlField()
        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormInf(explicit_control_field), 0.0000, 4)

    def test_GetPhysicalField(self):
        implicit_shape_field = self.implicit_shape_control.GePhysicalField()
        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormL2(implicit_shape_field), 20.976176963410882, 10)
        explicit_shape_field = self.explicit_shape_control.GePhysicalField()
        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormL2(explicit_shape_field), 20.976176963410882, 10)

    def test_MapGradient(self):
        physical_gradient = self.implicit_shape_control.GetEmptyField()
        constant_field_value = Kratos.Array3([1.0, 1.0, 1.0])
        Kratos.Expression.LiteralExpressionIO.SetData(physical_gradient, constant_field_value)
        self.explicit_shape_control.filter.GetIntegrationWeights(physical_gradient)
        implicit_mapped_gradient = self.implicit_shape_control.MapGradient({KratosOA.SHAPE: physical_gradient})
        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormL2(implicit_mapped_gradient), 27.84296340221239, 10)
        explicit_mapped_gradient = self.explicit_shape_control.MapGradient({KratosOA.SHAPE: physical_gradient})
        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormL2(explicit_mapped_gradient), 28.102970155390018, 10)

    def test_Update(self):
        update_field = self.implicit_shape_control.GetEmptyField()
        constant_field_value = Kratos.Array3([0.1, 0.1, 0.1])
        Kratos.Expression.LiteralExpressionIO.SetData(update_field, constant_field_value)
        self.implicit_shape_control.Update(update_field)
        implicit_control_field = self.implicit_shape_control.GetControlField()
        implicit_shape_field = self.implicit_shape_control.GePhysicalField()
        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormL2(implicit_control_field), 3.633180424916991, 10)
        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormL2(implicit_shape_field), 10.365298105786017, 10)

        update_field = self.explicit_shape_control.GetEmptyField()
        constant_field_value = Kratos.Array3([0.1, 0.1, 0.1])
        Kratos.Expression.LiteralExpressionIO.SetData(update_field, constant_field_value)
        self.explicit_shape_control.Update(update_field)
        explicit_control_field = self.explicit_shape_control.GetControlField()
        explicit_shape_field = self.explicit_shape_control.GePhysicalField()
        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormL2(explicit_control_field), 3.633180424916991, 10)
        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormL2(explicit_shape_field), 10.365298105786017, 10)

if __name__ == "__main__":
    kratos_unittest.main()
