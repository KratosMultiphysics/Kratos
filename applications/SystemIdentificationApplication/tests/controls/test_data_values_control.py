import KratosMultiphysics as Kratos
import KratosMultiphysics.StructuralMechanicsApplication as KratosSM
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.SystemIdentificationApplication.controls.data_values_control import DataValuesControl

class TestDataValuesControl_nodal_historical(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("Structure")
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3
        cls.model_part.AddNodalSolutionStepVariable(Kratos.TEMPERATURE)

        parameters = Kratos.Parameters("""{
            "control_variable_name"             : "TEMPERATURE",
            "control_variable_bounds"           : [-5.0, 5.0],
            "output_all_fields"                 : false,
            "container_type"                    : "nodal_historical",
            "filter_settings"                   : {
                    "filter_type": "explicit_filter",
                    "filter_function_type": "linear",
                    "max_nodes_in_filter_radius": 100000,
                    "echo_level": 4,
                    "filter_radius_settings": {
                        "filter_radius_type": "constant",
                        "filter_radius": 1e-10
                    }
            },
            "model_part_names": [
                {
                    "primal_model_part_name" : "Structure",
                    "adjoint_model_part_name": "Structure"
                }
            ]
        }""")

        cls.optimization_problem = OptimizationProblem()
        cls.temperature_control = DataValuesControl("test", cls.model, parameters, cls.optimization_problem)
        cls.optimization_problem.AddComponent(cls.temperature_control)

        with kratos_unittest.WorkFolderScope(".", __file__):
            Kratos.ModelPartIO("../responses/auxiliary_files_3/Structure", Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(cls.model_part)
        material_params = Kratos.Parameters("""{
            "properties": [
                {
                    "model_part_name": "Structure.Parts_Solid_full",
                    "properties_id": 1,
                    "Material": {
                        "constitutive_law" : {
                            "name" : "ConstitutiveLawsApplication.ThermalLinearPlaneStress"
                            },
                        "Variables"        : {
                            "THICKNESS"     : 0.02,
                            "DENSITY"       : 7850.0,
                            "YOUNG_MODULUS" : 2e6,
                            "POISSON_RATIO" : 0.29,
                            "THERMAL_EXPANSION_COEFFICIENT": 1e-5,
                            "REFERENCE_TEMPERATURE": 0.0
                        },
                        "Tables"           : null
                    }
                }
            ]
        }""")
        Kratos.ReadMaterialsUtility(cls.model).ReadMaterials(material_params)

        for node in cls.model_part.Nodes:
            node: Kratos.Node
            node.SetSolutionStepValue(Kratos.TEMPERATURE, -2.5)

        cls.temperature_control.Initialize()
        cls.initial_temperature = Kratos.Expression.NodalExpression(cls.model_part)
        Kratos.Expression.VariableExpressionIO.Read(cls.initial_temperature, Kratos.TEMPERATURE, 1)

    def setUp(self) -> None:
        Kratos.Expression.VariableExpressionIO.Write(self.initial_temperature, Kratos.TEMPERATURE, 1)

    @classmethod
    def tearDownClass(cls):
        with kratos_unittest.WorkFolderScope(".", __file__):
            DeleteFileIfExisting("Structure.time")

    def test_GetPhysicalKratosVariables(self):
        control_var = self.temperature_control.GetPhysicalKratosVariables()
        self.assertEqual(control_var, [Kratos.KratosGlobals.GetVariable("TEMPERATURE")])

    def test_GetControlField(self):
        control_field = self.temperature_control.GetControlField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormInf(control_field), 0.0, 4)

    def test_GetPhysicalField(self):
        temperature_field = self.temperature_control.GetPhysicalField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(temperature_field), 7.5, 4)

    def test_MapGradient(self):
        physical_gradient = self.temperature_control.GetEmptyField()
        for node in physical_gradient.GetContainer():
            node.SetValue(KratosSM.TEMPERATURE_SENSITIVITY, 1)
        Kratos.Expression.VariableExpressionIO.Read(physical_gradient, KratosSM.TEMPERATURE_SENSITIVITY, 0)
        mapped_gradient = self.temperature_control.MapGradient({Kratos.TEMPERATURE: physical_gradient})
        # physical = -2.5
        # ProjectBackward: phi = 0.5 - sin(asin(1-2*(physical - min)/delta)/3) = 0.3263518223
        # CalculateForwardProjectionGradient: d_physical/d_phi = (6*phi - 6*phi²)*delta = 13.19077862357725
        # d_J/d_physical = 1 (input given above)
        # BackwardFilterIntegratedField: d_J/d_physical * d_physical/d_phi -> d_J/d_control (mapped gradient)
        # For Integrated type: domain_size (of node_1) = element_area / num_nodes = 0.125 / 3 = 0.0416667
        # BackwardFilterIntegratedField: (1 * 13.19077862357725) / 0.0416667 = 316.57868697 (mapped gradient)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormInf(mapped_gradient), 316.57868697, 4)

    def test_Update(self):
        update_field = self.temperature_control.GetEmptyField()
        Kratos.Expression.LiteralExpressionIO.SetData(update_field, 0.25)
        temperature_field = self.temperature_control.GetPhysicalField()
        control_field = self.temperature_control.GetControlField()
        self.temperature_control.Update(update_field)
        control_field = self.temperature_control.GetControlField()
        temperature_field = self.temperature_control.GetPhysicalField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormInf(control_field), 0.25, 4)
        # ForwardFilter: control_update -> phi_update (Here, filter radius ~ 0. Therefore, both are the same = 0.25)
        # physical = -2.5
        # phi = 0.5 - sin(asin(1-2*(physical - min)/delta)/3) = 0.3263518223
        # phi_updated = phi_current + phi_update = 0.3263518223 + 0.25 = 0.5763518223
        # ProjectForward: phi_updated -> physical_updated
        # physical_updated = physical_min + phi_updated²*(3 - 2*phi_updated)* delta = 1.136375322
        self.assertAlmostEqual(Kratos.Expression.Utils.NormInf(temperature_field), 1.136375322, 6)

class TestDataValuesControl_condition(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("Structure")
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

        parameters = Kratos.Parameters("""{
            "control_variable_name"             : "PRESSURE",
            "control_variable_bounds"           : [-5.0, 5.0],
            "output_all_fields"                 : false,
            "container_type"                    : "condition",
            "filter_settings"                   : {
                    "filter_type": "explicit_filter",
                    "filter_function_type": "linear",
                    "max_nodes_in_filter_radius": 100000,
                    "echo_level": 4,
                    "filter_radius_settings": {
                        "filter_radius_type": "constant",
                        "filter_radius": 1e-10
                    }
            },
            "model_part_names": [
                {
                    "primal_model_part_name" : "Structure",
                    "adjoint_model_part_name": "Structure"
                }
            ]
        }""")

        cls.optimization_problem = OptimizationProblem()
        cls.pressure_control = DataValuesControl("test", cls.model, parameters, cls.optimization_problem)
        cls.optimization_problem.AddComponent(cls.pressure_control)

        with kratos_unittest.WorkFolderScope(".", __file__):
            Kratos.ModelPartIO("../responses/auxiliary_files_5/Structure", Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(cls.model_part)
        material_params = Kratos.Parameters("""{
            "properties": [
            {
                "model_part_name": "Structure.Parts_Solid",
                "properties_id": 1,
                "Material": {
                    "constitutive_law": {
                        "name": "LinearElastic3DLaw"
                    },
                    "Variables": {
                        "DENSITY": 7850.0,
                        "YOUNG_MODULUS": 20.0,
                        "POISSON_RATIO": 0.29
                    },
                    "Tables": null
                }
            }
    ]
        }""")
        Kratos.ReadMaterialsUtility(cls.model).ReadMaterials(material_params)

        for condition in cls.model_part.Conditions:
            condition: Kratos.Condition
            condition.SetValue(Kratos.PRESSURE, -2.5)

        cls.pressure_control.Initialize()
        cls.initial_pressure = Kratos.Expression.ConditionExpression(cls.model_part)
        Kratos.Expression.VariableExpressionIO.Read(cls.initial_pressure, Kratos.PRESSURE)

    def setUp(self) -> None:
        Kratos.Expression.VariableExpressionIO.Write(self.initial_pressure, Kratos.PRESSURE)

    @classmethod
    def tearDownClass(cls):
        with kratos_unittest.WorkFolderScope(".", __file__):
            DeleteFileIfExisting("Structure.time")

    def test_GetPhysicalKratosVariables(self):
        control_var = self.pressure_control.GetPhysicalKratosVariables()
        self.assertEqual(control_var, [Kratos.KratosGlobals.GetVariable("PRESSURE")])

    def test_GetControlField(self):
        control_field = self.pressure_control.GetControlField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormInf(control_field), 0.0, 4)

    def test_GetPhysicalField(self):
        pressure_field = self.pressure_control.GetPhysicalField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(pressure_field), 3.535533906, 4)

    def test_MapGradient(self):
        physical_gradient = self.pressure_control.GetEmptyField()
        for condition in physical_gradient.GetContainer():
            condition.SetValue(KratosSM.PRESSURE_SENSITIVITY, 1)
        Kratos.Expression.VariableExpressionIO.Read(physical_gradient, KratosSM.PRESSURE_SENSITIVITY)
        mapped_gradient = self.pressure_control.MapGradient({Kratos.PRESSURE: physical_gradient})
        # physical = -2.5
        # ProjectBackward: phi = 0.5 - sin(asin(1-2*(physical - min)/delta)/3) = 0.3263518223
        # CalculateForwardProjectionGradient: d_physical/d_phi = (6*phi - 6*phi²)*delta = 13.19077862357725
        # d_J/d_physical = 1 (input given above)
        # BackwardFilterIntegratedField: d_J/d_physical * d_physical/d_phi -> d_J/d_control (mapped gradient)
        # For Integrated type: domain_size (of condition_1) = condition_area = 0.25
        # BackwardFilterIntegratedField: (1 * 13.19077862357725) / 0.25 = 52.76311449 (mapped gradient)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormInf(mapped_gradient), 52.76311449, 4)

    def test_Update(self):
        update_field = self.pressure_control.GetEmptyField()
        Kratos.Expression.LiteralExpressionIO.SetData(update_field, 0.25)
        pressure_field = self.pressure_control.GetPhysicalField()
        control_field = self.pressure_control.GetControlField()
        self.pressure_control.Update(update_field)
        control_field = self.pressure_control.GetControlField()
        pressure_field = self.pressure_control.GetPhysicalField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormInf(control_field), 0.25, 4)
        # ForwardFilter: control_update -> phi_update (Here, filter radius ~ 0. Therefore, both are the same = 0.25)
        # physical = -2.5
        # phi = 0.5 - sin(asin(1-2*(physical - min)/delta)/3) = 0.3263518223
        # phi_updated = phi_current + phi_update = 0.3263518223 + 0.25 = 0.5763518223
        # ProjectForward: phi_updated -> physical_updated
        # physical_updated = physical_min + phi_updated²*(3 - 2*phi_updated)* delta = 1.136375322
        self.assertAlmostEqual(Kratos.Expression.Utils.NormInf(pressure_field), 1.136375322, 6)

if __name__ == "__main__":
    kratos_unittest.main()