
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.StructuralMechanicsApplication

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.controls.control_wrapper import ControlWrapper
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import ContainerEnum

class TestPropertiesControl(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.optimization_info = OptimizationInfo()
        cls.model_part = cls.model.CreateModelPart("Structure")
        Kratos.ModelPartIO("linear_element_test/Structure", Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(cls.model_part)

        material_settings = Kratos.Parameters("""{"Parameters": {"materials_filename": "linear_element_test/StructuralMaterials.json"}} """)
        Kratos.ReadMaterialsUtility(material_settings, cls.model)

        parameters = Kratos.Parameters("""{
            "name"                 : "density",
            "type"                 : "PropertiesControl",
            "settings"             : {
                "model_part_name"       : "Structure.structure",
                "control_variable_name" : "DENSITY"
            }
        }""")

        cls.properties_control_wrapper = ControlWrapper(cls.model, parameters, cls.optimization_info)
        cls.optimization_info.AddRoutine(cls.properties_control_wrapper)

    @classmethod
    def tearDownClass(cls):
        with kratos_unittest.WorkFolderScope("linear_element_test", __file__):
            DeleteFileIfExisting("Structure.time")

    def test_PropertiesControlInitialize(self):
        # running it twice to check whether the it only does the creation of specific properties once.
        self.optimization_info.Initialize()

        self.assertEqual(self.optimization_info["model_parts_with_element_specific_properties"], ["Structure.structure"])

        for element_i in self.model_part.Elements:
            for element_j in self.model_part.Elements:
                if element_i.Id != element_j.Id:
                    self.assertNotEqual(element_i.Properties, element_j.Properties)

    def test_PropertiesControl(self):
        self.optimization_info.Initialize()

        update_vector = Kratos.Vector()
        KratosOA.OptimizationUtils.GetContainerPropertiesVariableToVector(self.properties_control_wrapper.GetControl().GetModelPart().Elements, Kratos.DENSITY, update_vector)

        # run for 3 iterations
        for i in range(3):
            values = Kratos.Vector()
            KratosOA.OptimizationUtils.GetContainerPropertiesVariableToVector(self.properties_control_wrapper.GetControl().GetModelPart().Elements, Kratos.DENSITY, values)

            self.optimization_info.InitializeSolutionStep()

            if self.optimization_info["step"] > 1:
                for j, element in enumerate(self.properties_control_wrapper.GetControl().GetModelPart().Elements):
                    self.assertEqual(element.Properties[Kratos.DENSITY], values[j] + 1)

            self.properties_control_wrapper.GetControl().SetControlUpdatesVector(update_vector)

            self.optimization_info.FinalizeSolutionStep()

        self.optimization_info.Finalize()

    def test_GetModelPart(self):
        self.assertEqual(self.model_part.GetSubModelPart("structure"), self.properties_control_wrapper.GetControl().GetModelPart())

    def test_GetContainerType(self):
        self.assertEqual(self.properties_control_wrapper.GetControl().GetContainerType(), ContainerEnum.ELEMENT_PROPERTIES)

    def test_GetControlSensitivityVariable(self):
        self.assertEqual(self.properties_control_wrapper.GetControl().GetControlSensitivityVariable(), KratosOA.DENSITY_SENSITIVITY)

    def test_GetControlUpdateVariable(self):
        self.assertEqual(self.properties_control_wrapper.GetControl().GetControlUpdateVariable(), KratosOA.SCALAR_CONTROL_UPDATE)

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()