
import KratosMultiphysics as Kratos
from KratosMultiphysics.testing.utilities import ReadModelPart
import KratosMultiphysics.OptimizationApplication as KratosOA

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.controls.shape.vertex_morphing_shape_control import VertexMorphingShapeControl

class TestVMShapeControl(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.shell_model_part = cls.model.CreateModelPart("shell")
        cls.shell_model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

        with kratos_unittest.WorkFolderScope(".", __file__):
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
            ReadModelPart("shell",cls.shell_model_part)
            cls.implicit_shape_control.Initialize()
            cls.explicit_shape_control.Initialize()

            cls.solid_explicit = cls.model.CreateModelPart("solid_explicit")
            cls.solid_explicit.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

            solid_parameters = Kratos.Parameters("""{
                "controlled_model_part_names": ["solid_explicit.design"],
                "filter_settings": {
                    "type" : "surface_explicit_with_mesh_motion",
                    "radius": 0.5
                },
                "fixed_model_parts" : {
                    "solid_explicit.fixed" : [true,true,true]
                }
            }""")
            cls.explicit_solid_shape_control = VertexMorphingShapeControl("solid_explicit", cls.model, solid_parameters, cls.optimization_problem)
            cls.optimization_problem.AddComponent(cls.explicit_solid_shape_control)
            ReadModelPart("solid",cls.solid_explicit)
            cls.explicit_solid_shape_control.Initialize()

            cls.solid_implicit = cls.model.CreateModelPart("solid_implicit")
            cls.solid_implicit.ProcessInfo[Kratos.DOMAIN_SIZE] = 3
            solid_parameters = Kratos.Parameters("""{
                "controlled_model_part_names": ["solid_implicit.design"],
                "filter_settings": {
                    "type" : "surface_implicit_with_mesh_motion",
                    "radius": 0.5
                },
                "fixed_model_parts" : {
                    "solid_implicit.fixed" : [true,true,true]
                }
            }""")

            cls.implicit_solid_shape_control = VertexMorphingShapeControl("solid_implicit", cls.model, solid_parameters, cls.optimization_problem)
            cls.optimization_problem.AddComponent(cls.implicit_solid_shape_control)
            ReadModelPart("solid",cls.solid_implicit)
            cls.implicit_solid_shape_control.Initialize()


            cls.solid_bulk_surf_implicit = cls.model.CreateModelPart("solid_bulk_surf_implicit")
            cls.solid_bulk_surf_implicit.ProcessInfo[Kratos.DOMAIN_SIZE] = 3
            solid_parameters = Kratos.Parameters("""{
                "controlled_model_part_names": ["solid_bulk_surf_implicit.design"],
                "filter_settings": {
                    "type" : "bulk_surface_implicit",
                    "radius": 0.5
                },
                "fixed_model_parts" : {
                    "solid_bulk_surf_implicit.fixed" : [true,true,true]
                }
            }""")

            cls.implicit_bulk_surf_shape_control = VertexMorphingShapeControl("solid_bulk_surf_implicit", cls.model, solid_parameters, cls.optimization_problem)
            cls.optimization_problem.AddComponent(cls.implicit_bulk_surf_shape_control)
            ReadModelPart("solid",cls.solid_bulk_surf_implicit)
            cls.implicit_bulk_surf_shape_control.Initialize()

    @classmethod
    def tearDownClass(cls):
        with kratos_unittest.WorkFolderScope(".", __file__):
            DeleteFileIfExisting("shell.time")

    def test_GetControlField(self):
        implicit_control_field = self.implicit_shape_control.GetControlField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(implicit_control_field), 22.4268730405143, 10)
        explicit_control_field = self.explicit_shape_control.GetControlField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(explicit_control_field), 0.0000, 4)
        solid_explicit_control_field = self.explicit_solid_shape_control.GetControlField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(solid_explicit_control_field), 0.0000, 4)
        solid_implicit_control_field = self.implicit_solid_shape_control.GetControlField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(solid_implicit_control_field), 18.63464515358424, 4)
        solid_implicit_bulk_surf_control_field = self.implicit_bulk_surf_shape_control.GetControlField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(solid_implicit_bulk_surf_control_field), 213.27526425126743, 4)


    def test_GetPhysicalField(self):
        shape_field = self.implicit_shape_control.GetPhysicalField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(shape_field), 20.976176963410882, 10)
        shape_field = self.explicit_shape_control.GetPhysicalField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(shape_field), 20.976176963410882, 10)
        shape_field = self.explicit_solid_shape_control.GetPhysicalField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(shape_field), 4.06201920231798, 10)
        shape_field = self.implicit_solid_shape_control.GetPhysicalField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(shape_field), 4.06201920231798, 10)
        shape_field = self.implicit_bulk_surf_shape_control.GetPhysicalField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(shape_field), 4.24264068711927, 10)

    def test_MapGradient(self):
        physical_gradient = self.implicit_shape_control.GetEmptyField()
        constant_field_value = Kratos.Array3([1.0, 1.0, 1.0])
        Kratos.Expression.LiteralExpressionIO.SetData(physical_gradient, constant_field_value)
        self.explicit_shape_control.filter.GetIntegrationWeights(physical_gradient)
        implicit_mapped_gradient = self.implicit_shape_control.MapGradient({KratosOA.SHAPE: physical_gradient})
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(implicit_mapped_gradient), 27.84296340221239, 10)
        explicit_mapped_gradient = self.explicit_shape_control.MapGradient({KratosOA.SHAPE: physical_gradient})
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(explicit_mapped_gradient), 28.102970155390018, 10)

        solid_physical_gradient = self.explicit_solid_shape_control.GetEmptyField()
        Kratos.Expression.LiteralExpressionIO.SetData(solid_physical_gradient, constant_field_value)
        mapped_gradient = self.explicit_solid_shape_control.MapGradient({KratosOA.SHAPE: solid_physical_gradient})
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(mapped_gradient), 13.52774925846869, 10)

        solid_physical_gradient = self.implicit_solid_shape_control.GetEmptyField()
        Kratos.Expression.LiteralExpressionIO.SetData(solid_physical_gradient, constant_field_value)
        mapped_gradient = self.implicit_solid_shape_control.MapGradient({KratosOA.SHAPE: solid_physical_gradient})
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(mapped_gradient), 9.792958372369323, 10)

        solid_physical_gradient = self.implicit_bulk_surf_shape_control.GetEmptyField()
        Kratos.Expression.LiteralExpressionIO.SetData(solid_physical_gradient, constant_field_value)
        mapped_gradient = self.implicit_bulk_surf_shape_control.MapGradient({KratosOA.SHAPE: solid_physical_gradient})
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(mapped_gradient), 20.108624626448545, 10)

    def test_Update(self):
        update_field = self.implicit_shape_control.GetEmptyField()
        constant_field_value = Kratos.Array3([0.1, 0.1, 0.1])
        Kratos.Expression.LiteralExpressionIO.SetData(update_field, constant_field_value)
        self.implicit_shape_control.Update(update_field)
        implicit_control_field = self.implicit_shape_control.GetControlField()
        implicit_shape_field = self.implicit_shape_control.GetPhysicalField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(implicit_control_field), 3.633180424916991, 10)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(implicit_shape_field), 10.365298105786017, 10)

        update_field = self.explicit_shape_control.GetEmptyField()
        Kratos.Expression.LiteralExpressionIO.SetData(update_field, constant_field_value)
        self.explicit_shape_control.Update(update_field)
        explicit_control_field = self.explicit_shape_control.GetControlField()
        explicit_shape_field = self.explicit_shape_control.GetPhysicalField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(explicit_control_field), 3.633180424916991, 10)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(explicit_shape_field), 10.365298105786017, 10)

        update_field = self.explicit_solid_shape_control.GetEmptyField()
        Kratos.Expression.LiteralExpressionIO.SetData(update_field, constant_field_value)
        self.explicit_solid_shape_control.Update(update_field)
        control_field = self.explicit_solid_shape_control.GetControlField()
        shape_field = self.explicit_solid_shape_control.GetPhysicalField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(control_field), 0.6244997998398398, 10)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(shape_field), 4.355456348076514, 10)

        update_field = self.implicit_solid_shape_control.GetEmptyField()
        Kratos.Expression.LiteralExpressionIO.SetData(update_field, constant_field_value)
        self.implicit_solid_shape_control.Update(update_field)
        control_field = self.implicit_solid_shape_control.GetControlField()
        shape_field = self.implicit_solid_shape_control.GetPhysicalField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(control_field), 0.6244997998398398, 10)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(shape_field), 3.1351923726934183, 10)

        update_field = self.implicit_bulk_surf_shape_control.GetEmptyField()
        Kratos.Expression.LiteralExpressionIO.SetData(update_field, constant_field_value)
        self.implicit_bulk_surf_shape_control.Update(update_field)
        control_field = self.implicit_bulk_surf_shape_control.GetControlField()
        shape_field = self.implicit_bulk_surf_shape_control.GetPhysicalField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(control_field), 0.6480740698407862, 10)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(shape_field), 4.370590586296503, 10)

if __name__ == "__main__":
    kratos_unittest.main()
