
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.StructuralMechanicsApplication

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.controls.material_properties_control import MaterialPropertiesControl

class TestMaterialPropertiesControl(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.optimization_info = OptimizationInfo()
        cls.optimization_info.SetBufferSize(1)
        cls.optimization_info["step"] = 0
        cls.model_part = cls.model.CreateModelPart("Structure")
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3
        Kratos.ModelPartIO("linear_element_test/Structure", Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(cls.model_part)

        material_settings = Kratos.Parameters("""{"Parameters": {"materials_filename": "linear_element_test/StructuralMaterials.json"}} """)
        Kratos.ReadMaterialsUtility(material_settings, cls.model)

        parameters = Kratos.Parameters("""{
            "model_part_names"      : ["Structure.structure"],
            "control_variable_name" : "DENSITY"
        }""")

        cls.properties_control = MaterialPropertiesControl(cls.model, parameters, cls.optimization_info)

    @classmethod
    def tearDownClass(cls):
        with kratos_unittest.WorkFolderScope("linear_element_test", __file__):
            DeleteFileIfExisting("Structure.time")

    def setUp(self):
        self.optimization_info["step"] = 0

    def test_PropertiesControlInitialize(self):
        # running it twice to check whether the it only does the creation of specific properties once.
        self.properties_control.ExecuteInitialize()

        self.assertEqual(self.optimization_info["model_parts_with_element_specific_properties"], ["Structure.structure.Elements"])

        for element_i in self.model_part.Elements:
            for element_j in self.model_part.Elements:
                if element_i.Id != element_j.Id:
                    self.assertNotEqual(element_i.Properties, element_j.Properties)

    def test_PropertiesControl(self):
        self.properties_control.ExecuteInitialize()

        model_part = self.model_part.GetSubModelPart("structure")

        for element in model_part.Elements:
            element.Properties[Kratos.DENSITY] = element.Id

        update_vector = KratosOA.ElementPropertiesContainerVariableDataHolder(model_part)
        update_vector.ReadDataFromContainerVariable(Kratos.DENSITY)

        collective_update = KratosOA.CollectiveVariableDataHolder([update_vector])

        # run for 3 iterations
        for i in range(1, 4, 1):
            self.optimization_info.AdvanceSolutionStep()
            self.optimization_info["step"] = i
            self.properties_control.ExecuteInitializeSolutionStep()

            if i > 1:
                for element in model_part.Elements:
                    self.assertEqual(element.Properties[Kratos.DENSITY], element.Id * i)

            self.properties_control.UpdateControl(collective_update.Clone())

            self.properties_control.ExecuteFinalizeSolutionStep()

        self.properties_control.ExecuteFinalize()

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.TESTS_OUTPUTS)  # TESTS_OUTPUTS
    kratos_unittest.main()