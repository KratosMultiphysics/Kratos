import numpy as np
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
                    "max_items_in_bucket": 10,
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
        cls.initial_temperature = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(cls.model_part.Nodes, Kratos.TEMPERATURE)
        cls.initial_temperature.CollectData()

    def setUp(self) -> None:
        self.initial_temperature.CollectData()

    @classmethod
    def tearDownClass(cls):
        with kratos_unittest.WorkFolderScope(".", __file__):
            DeleteFileIfExisting("Structure.time")
            DeleteFileIfExisting("initial_nodal_historical_control_values.csv")
            DeleteFileIfExisting("../initial_nodal_historical_control_values.csv")

    def test_GetPhysicalKratosVariables(self):
        control_var = self.temperature_control.GetPhysicalKratosVariables()
        self.assertEqual(control_var, [Kratos.KratosGlobals.GetVariable("TEMPERATURE")])

    def test_GetControlField(self):
        control_field = self.temperature_control.GetControlField()
        self.assertAlmostEqual(np.linalg.norm(control_field.data), 0.0, 4)

    def test_GetPhysicalField(self):
        temperature_field = self.temperature_control.GetPhysicalField()
        self.assertAlmostEqual(np.linalg.norm(temperature_field.data), 7.5, 4)

    def test_MapGradient(self):
        physical_gradient = self.temperature_control.GetEmptyField()
        for node in physical_gradient.GetContainer():
            node.SetValue(KratosSM.TEMPERATURE_SENSITIVITY, 1)
        Kratos.TensorAdaptors.VariableTensorAdaptor(physical_gradient, KratosSM.TEMPERATURE_SENSITIVITY, copy=False).CollectData()

        mapped_gradient = self.temperature_control.MapGradient({Kratos.TEMPERATURE: physical_gradient})
        # physical = -2.5
        # ProjectBackward: phi = 0.5 - sin(asin(1-2*(physical - min)/delta)/3) = 0.3263518223
        # CalculateForwardProjectionGradient: d_physical/d_phi = (6*phi - 6*phi²)*delta = 13.19077862357725
        # d_J/d_physical = 1 (input given above)
        # BackwardFilterIntegratedField: d_J/d_physical * d_physical/d_phi -> d_J/d_control (mapped gradient)
        # For Integrated type: domain_size (of node_1) = element_area / num_nodes = 0.125 / 3 = 0.0416667
        # BackwardFilterIntegratedField: (1 * 13.19077862357725) / 0.0416667 = 316.57868697 (mapped gradient)
        self.assertAlmostEqual(np.linalg.norm(mapped_gradient.data, ord=np.inf), 316.57868697, 4)

    def test_Update(self):
        update_field = self.temperature_control.GetEmptyField()
        update_field.data[:] = 0.25
        temperature_field = self.temperature_control.GetPhysicalField()
        control_field = self.temperature_control.GetControlField()
        self.temperature_control.Update(update_field)
        control_field = self.temperature_control.GetControlField()
        temperature_field = self.temperature_control.GetPhysicalField()
        self.assertAlmostEqual(np.linalg.norm(control_field.data, ord=np.inf), 0.25, 4)
        # ForwardFilter: control_update -> phi_update (Here, filter radius ~ 0. Therefore, both are the same = 0.25)
        # physical = -2.5
        # IntervalBounder: min = -5, max = 5, delta = 10: GetBoundedTensorAdaptor: phi_b = (physical - min)/delta = 0.25
        # phi = 0.25 - sin(asin(1-2*(phi_b - min_clamper)/delta_clamper)/3) = 0.3263518223
        # phi_updated = phi_current + phi_update = 0.3263518223 + 0.25 = 0.5763518223
        # ProjectForward: phi_updated -> physical_updated
        # Clamper: physical_updated = min + phi_updated²*(3 - 2*phi_updated)* delta = 0.6136375322
        # GetUnboundedTensorAdaptor: physical_updated_unbounded = physical_min + physical_updated * delta = -5 + 0.6136375322 * 10 = 1.13637532
        self.assertAlmostEqual(np.linalg.norm(temperature_field.data, ord=np.inf), 1.13637532, 6)

class TestDataValuesControl_nodal_historical_with_initialization(TestDataValuesControl_nodal_historical):
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
                    "max_items_in_bucket": 10,
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
            ],
            "control_initialization_settings": {
                "initialize_from_file" : true,
                "filetype" : "csv",
                "file_name" : "initial_nodal_historical_control_values.csv"
            }
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
            node.SetSolutionStepValue(Kratos.TEMPERATURE, 5.0)

        # Create a CSV file with initial control values
        initial_control_values = np.ones(len(cls.model_part.Nodes))
        initial_control_values[:] = -2.5  # Set all initial control values to -2.5
        with open("initial_nodal_historical_control_values.csv", "w") as f:
            f.write("NodeId,TEMPERATURE\n")
            for idx, node in enumerate(cls.model_part.Nodes):
                f.write(f"{node.Id},{initial_control_values[idx]}\n")

        cls.temperature_control.Initialize()
        cls.initial_temperature = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(cls.model_part.Nodes, Kratos.TEMPERATURE)
        cls.initial_temperature.CollectData()

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
                    "max_items_in_bucket": 10,
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
        cls.initial_pressure = Kratos.TensorAdaptors.VariableTensorAdaptor(cls.model_part.Conditions, Kratos.PRESSURE)
        cls.initial_pressure.CollectData()

    def setUp(self) -> None:
        self.initial_pressure.CollectData()

    @classmethod
    def tearDownClass(cls):
        with kratos_unittest.WorkFolderScope(".", __file__):
            DeleteFileIfExisting("Structure.time")
            DeleteFileIfExisting("initial_condition_control_values.csv")
            DeleteFileIfExisting("../initial_condition_control_values.csv")

    def test_GetPhysicalKratosVariables(self):
        control_var = self.pressure_control.GetPhysicalKratosVariables()
        self.assertEqual(control_var, [Kratos.KratosGlobals.GetVariable("PRESSURE")])

    def test_GetControlField(self):
        control_field = self.pressure_control.GetControlField()
        self.assertAlmostEqual(np.linalg.norm(control_field.data, ord=np.inf), 0.0, 4)

    def test_GetPhysicalField(self):
        pressure_field = self.pressure_control.GetPhysicalField()
        self.assertAlmostEqual(np.linalg.norm(pressure_field.data), 3.535533906, 4)

    def test_MapGradient(self):
        physical_gradient = self.pressure_control.GetEmptyField()
        for condition in physical_gradient.GetContainer():
            condition.SetValue(KratosSM.PRESSURE_SENSITIVITY, 1)
        Kratos.TensorAdaptors.VariableTensorAdaptor(physical_gradient, KratosSM.PRESSURE_SENSITIVITY, copy=False).CollectData()
        mapped_gradient = self.pressure_control.MapGradient({Kratos.PRESSURE: physical_gradient})
        # physical = -2.5
        # ProjectBackward: phi = 0.5 - sin(asin(1-2*(physical - min)/delta)/3) = 0.3263518223
        # CalculateForwardProjectionGradient: d_physical/d_phi = (6*phi - 6*phi²)*delta = 13.19077862357725
        # d_J/d_physical = 1 (input given above)
        # BackwardFilterIntegratedField: d_J/d_physical * d_physical/d_phi -> d_J/d_control (mapped gradient)
        # For Integrated type: domain_size (of condition_1) = condition_area = 0.25
        # BackwardFilterIntegratedField: (1 * 13.19077862357725) / 0.25 = 52.76311449 (mapped gradient)
        self.assertAlmostEqual(np.linalg.norm(mapped_gradient.data, ord=np.inf), 52.76311449, 4)

    def test_Update(self):
        update_field = self.pressure_control.GetEmptyField()
        update_field.data[:] = 0.25
        pressure_field = self.pressure_control.GetPhysicalField()
        control_field = self.pressure_control.GetControlField()
        self.pressure_control.Update(update_field)
        control_field = self.pressure_control.GetControlField()
        pressure_field = self.pressure_control.GetPhysicalField()
        self.assertAlmostEqual(np.linalg.norm(control_field.data, ord=np.inf), 0.25, 4)
        # ForwardFilter: control_update -> phi_update (Here, filter radius ~ 0. Therefore, both are the same = 0.25)
        # physical = -2.5
        # phi = 0.5 - sin(asin(1-2*(physical - min)/delta)/3) = 0.3263518223
        # phi_updated = phi_current + phi_update = 0.3263518223 + 0.25 = 0.5763518223
        # ProjectForward: phi_updated -> physical_updated
        # physical_updated = physical_min + phi_updated²*(3 - 2*phi_updated)* delta = 1.136375322
        self.assertAlmostEqual(np.linalg.norm(pressure_field.data, ord=np.inf), 1.136375322, 6)

class TestDataValuesControl_condition_with_initialization(TestDataValuesControl_condition):
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
                    "max_items_in_bucket": 10,
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
            ],
            "control_initialization_settings": {
                "initialize_from_file" : true,
                "filetype" : "csv",
                "file_name" : "initial_condition_control_values.csv"
            }
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
            condition.SetValue(Kratos.PRESSURE, 10.0)

        # Create a CSV file with initial control values
        initial_control_values = np.ones(len(cls.model_part.Conditions))
        initial_control_values[:] = -2.5  # Set all initial control values to -2.5
        with open("initial_condition_control_values.csv", "w") as f:
            f.write("ConditionId,PRESSURE\n")
            for idx, cond in enumerate(cls.model_part.Conditions):
                f.write(f"{cond.Id},{initial_control_values[idx]}\n")

        cls.pressure_control.Initialize()
        cls.initial_pressure = Kratos.TensorAdaptors.VariableTensorAdaptor(cls.model_part.Conditions, Kratos.PRESSURE)
        cls.initial_pressure.CollectData()

if __name__ == "__main__":
    kratos_unittest.main()