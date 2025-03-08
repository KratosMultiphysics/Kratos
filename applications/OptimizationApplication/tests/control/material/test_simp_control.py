
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.controls.material.simp_control import SimpControl

class TestSimpControl(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("shell")
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

        parameters = Kratos.Parameters("""{
            "controlled_model_part_names": [
                "shell"
            ],
            "output_all_fields": true,
            "list_of_materials": [
                {
                    "density": 0.0,
                    "young_modulus": 0.0
                },
                {
                    "density": 7850.0,
                    "young_modulus": 206900000000.0
                }
            ],
            "density_projection_settings": {
                "type": "adaptive_sigmoidal_projection",
                "initial_value": 2,
                "max_value": 50,
                "increase_fac": 1.05,
                "update_period": 25
            },
            "young_modulus_projection_settings": {
                "type": "adaptive_sigmoidal_projection",
                "initial_value": 2,
                "max_value": 50,
                "increase_fac": 1.05,
                "update_period": 25,
                "penalty_factor": 3
            },
            "filter_settings": {
                "filter_type": "explicit_filter",
                "filter_function_type": "linear",
                "max_items_in_bucket": 10,
                "echo_level": 0,
                "filter_radius_settings": {
                    "filter_radius_type": "constant",
                    "filter_radius": 0.2
                },
                "filtering_boundary_conditions": {
                    "damping_type": "nearest_entity",
                    "damping_function_type": "cosine",
                    "damped_model_part_settings": {}
                }
            }
        }""")

        cls.optimization_problem = OptimizationProblem()
        cls.simp_control = SimpControl("test", cls.model, parameters, cls.optimization_problem)
        cls.optimization_problem.AddComponent(cls.simp_control)

        with kratos_unittest.WorkFolderScope(".", __file__):
            Kratos.ModelPartIO("../thickness/Thick_2x2_Shell", Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(cls.model_part)
        material_params = Kratos.Parameters("""{
            "properties": [
                {
                    "model_part_name": "shell",
                    "properties_id": 1,
                    "Material": {
                        "Variables": {
                            "DENSITY"       : 3925.0,
                            "POISSON_RATIO" : 0.29
                        },
                        "Tables": {}
                    }
                }
            ]
        }""")
        Kratos.ReadMaterialsUtility(cls.model).ReadMaterials(material_params)

        cls.simp_control.Initialize()
        cls.initial_density = Kratos.Expression.ElementExpression(cls.model_part)
        KratosOA.PropertiesVariableExpressionIO.Read(cls.initial_density, Kratos.DENSITY)

    def setUp(self) -> None:
        KratosOA.PropertiesVariableExpressionIO.Write(self.initial_density, Kratos.DENSITY)

    @classmethod
    def tearDownClass(cls):
        with kratos_unittest.WorkFolderScope(".", __file__):
            DeleteFileIfExisting("../thickness/Thick_2x2_Shell.time")

    def test_GetControlField(self):
        control_field = self.simp_control.GetControlField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormInf(control_field), 0.0, 4)

    def test_MapGradient(self):
        physical_rho_gradient = self.simp_control.GetEmptyField()
        physical_E_gradient = self.simp_control.GetEmptyField()
        for element in physical_rho_gradient.GetContainer():
            element.SetValue(KratosOA.DENSITY_SENSITIVITY, element.GetGeometry().DomainSize())
            element.SetValue(KratosOA.YOUNG_MODULUS_SENSITIVITY, element.GetGeometry().DomainSize()/1e+10)
        Kratos.Expression.VariableExpressionIO.Read(physical_rho_gradient, KratosOA.DENSITY_SENSITIVITY)
        Kratos.Expression.VariableExpressionIO.Read(physical_E_gradient, KratosOA.YOUNG_MODULUS_SENSITIVITY)
        mapped_gradient = self.simp_control.MapGradient({Kratos.DENSITY: physical_rho_gradient, Kratos.YOUNG_MODULUS: physical_E_gradient})
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(mapped_gradient), 15731.035, 4)

    def test_Update(self):
        rho_field = Kratos.Expression.ElementExpression(self.model_part)
        KratosOA.PropertiesVariableExpressionIO.Read(rho_field, Kratos.DENSITY)
        E_field = Kratos.Expression.ElementExpression(self.model_part)
        KratosOA.PropertiesVariableExpressionIO.Read(E_field, Kratos.YOUNG_MODULUS)

        # update with empty field
        update_field = self.simp_control.GetEmptyField()
        self.simp_control.Update(update_field)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormInf(self.simp_control.GetControlField()), 0.0)
        updated_rho_field = Kratos.Expression.ElementExpression(self.model_part)
        KratosOA.PropertiesVariableExpressionIO.Read(updated_rho_field, Kratos.DENSITY)
        updated_E_field = Kratos.Expression.ElementExpression(self.model_part)
        KratosOA.PropertiesVariableExpressionIO.Read(updated_E_field, Kratos.YOUNG_MODULUS)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormInf(updated_rho_field - rho_field), 0.0)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormInf(updated_E_field - E_field), 0.0)

        Kratos.Expression.LiteralExpressionIO.SetData(update_field, 0.75)
        self.simp_control.Update(update_field)
        control_field = self.simp_control.GetControlField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormInf(control_field), 0.75, 4)

        KratosOA.PropertiesVariableExpressionIO.Read(updated_rho_field, Kratos.DENSITY)
        KratosOA.PropertiesVariableExpressionIO.Read(updated_E_field, Kratos.YOUNG_MODULUS)
        self.assertTrue(Kratos.Expression.Utils.NormInf(updated_rho_field - rho_field) > 0.0)
        self.assertTrue(Kratos.Expression.Utils.NormInf(updated_E_field - E_field) > 0.0)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(updated_rho_field), 14955.413791112203)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(updated_E_field), 357673554186.7966)

    def test_AdaptiveBeta(self):
        parameters = Kratos.Parameters("""{
            "controlled_model_part_names": ["shell"],
            "list_of_materials": [
                {
                    "density": 0.0,
                    "young_modulus": 0.0
                },
                {
                    "density": 7850.0,
                    "young_modulus": 206900000000.0
                }
            ],
            "filter_settings": {
                "filter_type": "explicit_filter",
                "filter_function_type": "linear",
                "max_items_in_bucket": 100000,
                "echo_level": 0,
                "filter_radius_settings": {
                    "filter_radius_type": "constant",
                    "filter_radius": 0.2
                },
                "filtering_boundary_conditions": {
                    "damping_type": "nearest_entity",
                    "damping_function_type": "cosine",
                    "damped_model_part_settings": {}
                }
            },
            "density_projection_settings": {
                "type": "adaptive_sigmoidal_projection",
                "initial_value": 0.01,
                "max_value"    : 30,
                "increase_fac" : 1.05,
                "update_period": 3
            },
            "young_modulus_projection_settings": {
                "type": "adaptive_sigmoidal_projection",
                "initial_value": 0.01,
                "max_value"    : 30,
                "increase_fac" : 1.05,
                "update_period": 3,
                "penalty_factor": 3
            }
        }""")
        simp_control = SimpControl("test_adaptive", self.model, parameters, self.optimization_problem)
        self.optimization_problem.AddComponent(simp_control)
        simp_control.Initialize()

        self.assertAlmostEqual(simp_control.density_projection.beta, 0.01)

        control_field = simp_control.GetControlField()
        simp_control.Update(control_field)
        for i in range(20):
            control_field *= 1.2
            simp_control.Update(control_field)
            self.optimization_problem.AdvanceStep()

        self.assertAlmostEqual(simp_control.density_projection.beta, 0.01 * (1.05 ** 7))

    def test_NonRemovalOfMaterials(self):
        parameters = Kratos.Parameters("""{
            "controlled_model_part_names": ["shell"],
            "list_of_materials": [
                {
                    "density": 3000,
                    "young_modulus": 106900000000
                },
                {
                    "density": 7850.0,
                    "young_modulus": 206900000000.0
                }
            ],
            "filter_settings": {
                "filter_type": "explicit_filter",
                "filter_function_type": "linear",
                "max_items_in_bucket": 100000,
                "echo_level": 0,
                "filter_radius_settings": {
                    "filter_radius_type": "constant",
                    "filter_radius": 0.2
                },
                "filtering_boundary_conditions": {
                    "damping_type": "nearest_entity",
                    "damping_function_type": "cosine",
                    "damped_model_part_settings": {}
                }
            },
            "density_projection_settings": {
                "type": "adaptive_sigmoidal_projection",
                "initial_value": 0.01,
                "max_value"    : 30,
                "increase_fac" : 1.05,
                "update_period": 3
            },
            "young_modulus_projection_settings": {
                "type": "adaptive_sigmoidal_projection",
                "initial_value": 0.01,
                "max_value"    : 30,
                "increase_fac" : 1.05,
                "update_period": 3,
                "penalty_factor": 3
            }
        }""")

        simp_control = SimpControl("test_non_removal", self.model, parameters, self.optimization_problem)
        self.optimization_problem.AddComponent(simp_control)

        simp_control.Initialize()

        control_field = simp_control.GetControlField()
        simp_control.Update(control_field)
        for i in range(20):
            control_field *= 1.2
            simp_control.Update(control_field)

        density = Kratos.Expression.ElementExpression(self.model_part)
        KratosOA.PropertiesVariableExpressionIO.Read(density, Kratos.DENSITY)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(density - 3925), 0.0)

if __name__ == "__main__":
    kratos_unittest.main()