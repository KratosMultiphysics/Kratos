
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.StructuralMechanicsApplication

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.controls.control_wrapper import ControlWrapper
from KratosMultiphysics.OptimizationApplication.utilities.container_data import ContainerData

class TestMaterialPropertiesControl(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.optimization_info = OptimizationInfo()
        cls.optimization_info.SetBufferSize(1)
        cls.optimization_info["step"] = 0
        cls.model_part = cls.model.CreateModelPart("Structure")
        Kratos.ModelPartIO("linear_element_test/Structure", Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(cls.model_part)

        material_settings = Kratos.Parameters("""{"Parameters": {"materials_filename": "linear_element_test/StructuralMaterials.json"}} """)
        Kratos.ReadMaterialsUtility(material_settings, cls.model)

        parameters = Kratos.Parameters("""{
            "name"                 : "density",
            "type"                 : "MaterialPropertiesControl",
            "settings"             : {
                "model_part_names"      : ["Structure.structure"],
                "control_variable_name" : "DENSITY"
            }
        }""")

        cls.properties_control_wrapper = ControlWrapper(cls.model, parameters, cls.optimization_info)
        cls.optimization_info.AddOptimizationRoutine(cls.properties_control_wrapper)

    @classmethod
    def tearDownClass(cls):
        with kratos_unittest.WorkFolderScope("linear_element_test", __file__):
            DeleteFileIfExisting("Structure.time")

    def setUp(self):
        self.optimization_info["step"] = 0

    def test_PropertiesControlInitialize(self):
        # running it twice to check whether the it only does the creation of specific properties once.
        self.properties_control_wrapper.Initialize()

        self.assertEqual(self.optimization_info["model_parts_with_element_specific_properties"], ["Structure.structure.Elements"])

        for element_i in self.model_part.Elements:
            for element_j in self.model_part.Elements:
                if element_i.Id != element_j.Id:
                    self.assertNotEqual(element_i.Properties, element_j.Properties)

    def test_PropertiesControl(self):
        self.properties_control_wrapper.Initialize()

        for element in self.model_part.GetSubModelPart("structure").Elements:
            element.Properties[Kratos.DENSITY] = element.Id

        update_vector = ContainerData(self.properties_control_wrapper.GetControl().GetModelParts()[0], self.properties_control_wrapper.GetControl().GetContainerType())
        update_vector.ReadDataFromContainer(Kratos.DENSITY)

        # run for 3 iterations
        for i in range(1, 4, 1):
            values = Kratos.Vector()
            KratosOA.OptimizationUtils.GetContainerPropertiesVariableToVector(self.properties_control_wrapper.GetControl().GetModelParts()[0].Elements, Kratos.DENSITY, values)

            self.optimization_info.AdvanceSolutionStep()
            self.optimization_info["step"] = i
            self.properties_control_wrapper.InitializeSolutionStep()

            if i > 1:
                for element in self.properties_control_wrapper.GetControl().GetModelParts()[0].Elements:
                    self.assertEqual(element.Properties[Kratos.DENSITY], element.Id * i)

            self.properties_control_wrapper.SetControlUpdate(update_vector.Clone())

            self.properties_control_wrapper.FinalizeSolutionStep()

        self.properties_control_wrapper.Finalize()

    def test_GetModelPart(self):
        self.assertEqual(self.model_part.GetSubModelPart("structure"), self.properties_control_wrapper.GetControl().GetModelParts()[0])

    def test_GetContainerType(self):
        self.assertEqual(self.properties_control_wrapper.GetControl().GetContainerType(), ContainerData.ContainerEnum.ELEMENT_PROPERTIES)

    def test_GetControlSensitivityVariable(self):
        self.assertEqual(self.properties_control_wrapper.GetControl().GetControlSensitivityVariable(), KratosOA.DENSITY_SENSITIVITY)

    def test_GetControlUpdateVariable(self):
        self.assertEqual(self.properties_control_wrapper.GetControl().GetControlUpdateVariable(), KratosOA.SCALAR_CONTROL_UPDATE)

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.TESTS_OUTPUTS)  # TESTS_OUTPUTS
    kratos_unittest.main()