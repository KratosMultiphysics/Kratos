import numpy as np
import KratosMultiphysics as Kratos
import KratosMultiphysics.StructuralMechanicsApplication as KratosSM
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.SystemIdentificationApplication.controls.material_properties_control import MaterialPropertiesControl

class TestMaterialPropertiesControl(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("Structure")
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

        parameters = Kratos.Parameters("""{
            "control_variable_name"             : "YOUNG_MODULUS",
            "control_variable_bounds"           : [1e6, 3e6],
            "output_all_fields"                 : false,
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
        cls.young_modulus_control = MaterialPropertiesControl("test", cls.model, parameters, cls.optimization_problem)
        cls.optimization_problem.AddComponent(cls.young_modulus_control)

        with kratos_unittest.WorkFolderScope(".", __file__):
            Kratos.ModelPartIO("../responses/auxiliary_files_3/Structure", Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(cls.model_part)
        material_params = Kratos.Parameters("""{
            "properties": [
                {
                    "model_part_name": "Structure.Parts_Solid_full",
                    "properties_id": 1,
                    "Material": {
                        "constitutive_law" : {
                            "name" : "LinearElasticPlaneStress2DLaw"
                            },
                        "Variables"        : {
                            "THICKNESS"     : 0.02,
                            "DENSITY"       : 7850.0,
                            "YOUNG_MODULUS" : 2e6,
                            "POISSON_RATIO" : 0.29
                        },
                        "Tables"           : null
                    }
                }
            ]
        }""")
        Kratos.ReadMaterialsUtility(cls.model).ReadMaterials(material_params)

        cls.young_modulus_control.Initialize()
        cls.initial_young_modulus = KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor(cls.model_part.Elements, Kratos.YOUNG_MODULUS)
        cls.initial_young_modulus.CollectData()

    def setUp(self) -> None:
        self.initial_young_modulus.CollectData()

    @classmethod
    def tearDownClass(cls):
        with kratos_unittest.WorkFolderScope(".", __file__):
            DeleteFileIfExisting("Structure.time")
            DeleteFileIfExisting("initial_material_properties_control_values.csv")
            DeleteFileIfExisting("../initial_material_properties_control_values.csv")

    def test_GetPhysicalKratosVariables(self):
        control_var = self.young_modulus_control.GetPhysicalKratosVariables()
        self.assertEqual(control_var, [Kratos.KratosGlobals.GetVariable("YOUNG_MODULUS")])

    def test_GetControlField(self):
        control_field = self.young_modulus_control.GetControlField()
        self.assertAlmostEqual(np.linalg.norm(control_field.data), 0.0, 4)

    def test_GetPhysicalField(self):
        young_modulus_field = self.young_modulus_control.GetPhysicalField()
        self.assertAlmostEqual(np.linalg.norm(young_modulus_field.data), 5656854.2494924, 4)

    def test_MapGradient(self):
        physical_gradient = self.young_modulus_control.GetEmptyField()
        for element in physical_gradient.GetContainer():
            element.SetValue(KratosSM.YOUNG_MODULUS_SENSITIVITY, 1)
        Kratos.TensorAdaptors.VariableTensorAdaptor(physical_gradient, KratosSM.YOUNG_MODULUS_SENSITIVITY, copy=False).CollectData()

        mapped_gradient = self.young_modulus_control.MapGradient({Kratos.YOUNG_MODULUS: physical_gradient})
        # physical = 2e6
        # interval bounder: min = 1e6, max = 3e6, delta = 2e6: GetBoundedTensorAdaptor: phi_b = (physical - min)/delta = 0.5
        # Clamper: ProjectBackward: phi = 0.5 - sin(asin(1-2*(phi_b - min_clamper)/delta_clamper)/3) = 0.5
        # Clamper: CalculateForwardProjectionGradient: d_physical/d_phi = (6*phi - 6*phi²)*delta_physical = 1.5 * 2e6 = 3e6
        # d_J/d_physical = 1 (input given above)
        # BackwardFilterIntegratedField: d_J/d_physical * d_physical/d_phi -> d_J/d_control (mapped gradient)
        # For Integrated type: domain_size  = element_area  = 0.125
        # BackwardFilterIntegratedField: (1 * 3e6) / 0.125 = 24000000.0 (mapped gradient)
        self.assertAlmostEqual(np.linalg.norm(mapped_gradient.data, ord=np.inf), 24000000.0, 4)

    def test_Update(self):
        update_field = self.young_modulus_control.GetEmptyField()
        update_field.data[:] = 0.25
        young_modulus_field = self.young_modulus_control.GetPhysicalField()
        control_field = self.young_modulus_control.GetControlField()
        self.young_modulus_control.Update(update_field)
        control_field = self.young_modulus_control.GetControlField()
        young_modulus_field = self.young_modulus_control.GetPhysicalField()
        self.assertAlmostEqual(np.linalg.norm(control_field.data, ord=np.inf), 0.25, 4)
        # ForwardFilter: control_update -> phi_update (Here, filter radius ~ 0. Therefore, both are the same = 0.25)
        # physical = 2e6
        # IntervalBounder: min = 1e6, max = 3e6, delta = 2e6: GetBoundedTensorAdaptor: phi_b = (physical - min)/delta = 0.5
        # phi = 0.5 - sin(asin(1-2*(phi_b - min_clamper)/delta_clamper)/3) = 0.5
        # phi_updated = phi_current + phi_update = 0.5 + 0.25 = 0.75
        # ProjectForward: phi_updated -> physical_updated
        # Clamper: physical_updated = min + phi_updated²*(3 - 2*phi_updated)* delta = 0.84375
        # GetUnboundedTensorAdaptor: physical_updated_unbounded = physical_min + physical_updated * delta = 1e6 + 0.84375 * 2e6 = 2.6875e6
        self.assertAlmostEqual(np.linalg.norm(young_modulus_field.data, ord=np.inf), 2.6875e6, 9)

class TestMaterialPropertiesControl_with_initialization(TestMaterialPropertiesControl):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("Structure")
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

        parameters = Kratos.Parameters("""{
            "control_variable_name"             : "YOUNG_MODULUS",
            "control_variable_bounds"           : [1e6, 3e6],
            "output_all_fields"                 : false,
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
                "file_name" : "initial_material_properties_control_values.csv"
            }
        }""")

        cls.optimization_problem = OptimizationProblem()
        cls.young_modulus_control = MaterialPropertiesControl("test", cls.model, parameters, cls.optimization_problem)
        cls.optimization_problem.AddComponent(cls.young_modulus_control)

        with kratos_unittest.WorkFolderScope(".", __file__):
            Kratos.ModelPartIO("../responses/auxiliary_files_3/Structure", Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(cls.model_part)
        material_params = Kratos.Parameters("""{
            "properties": [
                {
                    "model_part_name": "Structure.Parts_Solid_full",
                    "properties_id": 1,
                    "Material": {
                        "constitutive_law" : {
                            "name" : "LinearElasticPlaneStress2DLaw"
                            },
                        "Variables"        : {
                            "THICKNESS"     : 0.02,
                            "DENSITY"       : 7850.0,
                            "YOUNG_MODULUS" : 1.1e6,
                            "POISSON_RATIO" : 0.29
                        },
                        "Tables"           : null
                    }
                }
            ]
        }""")
        Kratos.ReadMaterialsUtility(cls.model).ReadMaterials(material_params)

        # Create a CSV file with initial control values
        initial_control_values = np.ones(len(cls.model_part.Elements))
        initial_control_values[:] = 2e6  # Set all initial control values to 2e6
        with open("initial_material_properties_control_values.csv", "w") as f:
            f.write("ElementId,YOUNG_MODULUS\n")
            for idx, element in enumerate(cls.model_part.Elements):
                f.write(f"{element.Id},{initial_control_values[idx]}\n")

        cls.young_modulus_control.Initialize()
        cls.initial_young_modulus = KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor(cls.model_part.Elements, Kratos.YOUNG_MODULUS)
        cls.initial_young_modulus.CollectData()

if __name__ == "__main__":
    kratos_unittest.main()