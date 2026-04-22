
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.StructuralMechanicsApplication

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.OptimizationApplication.controls.material.material_properties_control import MaterialPropertiesControl

class TestMaterialPropertiesControl(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("Structure")
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

        with kratos_unittest.WorkFolderScope("../../responses_tests/linear_strain_energy_test", __file__):
            Kratos.ModelPartIO("Structure", Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(cls.model_part)
            material_settings = Kratos.Parameters("""{"Parameters": {"materials_filename": "StructuralMaterials.json"}} """)
            Kratos.ReadMaterialsUtility(material_settings, cls.model)

        parameters = Kratos.Parameters("""{
            "model_part_names"      : ["Structure.structure"],
            "control_variable_name" : "DENSITY"
        }""")

        cls.properties_control = MaterialPropertiesControl("test", cls.model, parameters)

    @classmethod
    def tearDownClass(cls):
        with kratos_unittest.WorkFolderScope("../../responses_tests/linear_strain_energy_test", __file__):
            DeleteFileIfExisting("Structure.time")

    def test_PropertiesControlInitialize(self):
        # running it twice to check whether the it only does the creation of specific properties once.
        self.properties_control.Initialize()

        for element_i in self.model_part.Elements:
            for element_j in self.model_part.Elements:
                if element_i.Id != element_j.Id:
                    self.assertNotEqual(element_i.Properties, element_j.Properties)

    def test_PropertiesControl(self):
        self.properties_control.Initialize()

        model_part = self.model_part.GetSubModelPart("structure")

        for element in model_part.Elements:
            element.Properties[Kratos.DENSITY] = element.Id

        update_vector = self.properties_control.GetEmptyField()
        KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor(update_vector, Kratos.DENSITY, copy=False).CollectData()

        # run for 3 iterations
        for i in range(1, 4, 1):
            current_ta = Kratos.TensorAdaptors.DoubleTensorAdaptor(update_vector)
            current_ta.data *= i
            self.properties_control.Update(current_ta)

            for element in model_part.Elements:
                self.assertEqual(element.Properties[Kratos.DENSITY], element.Id * i)

        self.properties_control.Finalize()

    def test_PropertiesControlUpdate(self):
        self.properties_control.Initialize()

        control_model_part = self.model_part.GetSubModelPart("structure")

        with self.assertRaises(RuntimeError):
            self.properties_control.Update(Kratos.TensorAdaptors.VariableTensorAdaptor(control_model_part.Conditions, Kratos.DENSITY))

        temp = KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor(self.model_part.Elements, Kratos.DENSITY)
        temp.CollectData()

        with self.assertRaises(RuntimeError):
            self.properties_control.Update(temp)

        temp = KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor(control_model_part.Elements, Kratos.DENSITY)
        temp.CollectData()

        self.properties_control.Update(temp)

    def test_PropertiesControlMapGradient(self):
        self.properties_control.Initialize()

        control_model_part = self.model_part.GetSubModelPart("structure")

        with self.assertRaises(RuntimeError):
            self.properties_control.MapGradient({Kratos.DENSITY: Kratos.TensorAdaptors.VariableTensorAdaptor(control_model_part.Nodes, Kratos.DENSITY)})

        with self.assertRaises(RuntimeError):
            self.properties_control.MapGradient({Kratos.DENSITY: Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Elements, Kratos.DENSITY)})

        with self.assertRaises(RuntimeError):
            self.properties_control.MapGradient({Kratos.THICKNESS: Kratos.TensorAdaptors.VariableTensorAdaptor(control_model_part.Elements, Kratos.DENSITY)})

        with self.assertRaises(RuntimeError):
            self.properties_control.MapGradient({
                Kratos.DENSITY: Kratos.TensorAdaptors.VariableTensorAdaptor(control_model_part.Elements, Kratos.DENSITY),
                Kratos.THICKNESS: Kratos.TensorAdaptors.VariableTensorAdaptor(control_model_part.Elements, Kratos.DENSITY)})

        self.assertTrue(self.properties_control.MapGradient({Kratos.DENSITY: Kratos.TensorAdaptors.VariableTensorAdaptor(control_model_part.Elements, Kratos.DENSITY)}).GetContainer(), control_model_part.Elements)

    def test_GetControlFiield(self):
        self.properties_control.Initialize()
        field = self.properties_control.GetControlField()
        KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor(field, Kratos.YOUNG_MODULUS, copy=False).StoreData()
        for element in self.model_part.Elements:
            properties = element.Properties
            self.assertEqual(properties[Kratos.DENSITY], properties[Kratos.YOUNG_MODULUS])

if __name__ == "__main__":
    kratos_unittest.main()