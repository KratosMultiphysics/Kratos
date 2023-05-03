
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.StructuralMechanicsApplication

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.controls.material.material_properties_control import MaterialPropertiesControl

class TestMaterialPropertiesControl(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.optimization_problem = OptimizationProblem()
        cls.model_part = cls.model.CreateModelPart("Structure")
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3
        Kratos.ModelPartIO("linear_strain_energy_test/Structure", Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(cls.model_part)

        material_settings = Kratos.Parameters("""{"Parameters": {"materials_filename": "linear_strain_energy_test/StructuralMaterials.json"}} """)
        Kratos.ReadMaterialsUtility(material_settings, cls.model)

        parameters = Kratos.Parameters("""{
            "combined_output_model_part_name": "<CONTROL_NAME>",
            "model_part_names"      : ["Structure.structure"],
            "control_variable_name" : "DENSITY"
        }""")

        cls.properties_control = MaterialPropertiesControl(cls.model, parameters, cls.optimization_problem)
        cls.optimization_problem.AddControl("test", cls.properties_control)

    @classmethod
    def tearDownClass(cls):
        with kratos_unittest.WorkFolderScope("linear_strain_energy_test", __file__):
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

        update_vector = self.properties_control.GetEmptyControlField()
        update_vector.Read(Kratos.DENSITY)

        # run for 3 iterations
        for i in range(1, 4, 1):
            self.optimization_problem.AdvanceStep()

            self.properties_control.Update(update_vector.Clone() * i)

            for element in model_part.Elements:
                self.assertEqual(element.Properties[Kratos.DENSITY], element.Id * i)

        self.properties_control.Finalize()

    def test_PropertiesControlUpdate(self):
        self.properties_control.Initialize()

        control_model_part = self.model_part.GetSubModelPart("test")

        with self.assertRaises(RuntimeError):
            self.properties_control.Update(Kratos.ContainerExpression.ConditionNonHistoricalExpression(control_model_part))

        temp = KratosOA.ContainerExpression.ElementPropertiesExpression(self.model_part)
        temp.Read(Kratos.DENSITY)

        with self.assertRaises(RuntimeError):
            self.properties_control.Update(temp)

        temp = KratosOA.ContainerExpression.ElementPropertiesExpression(control_model_part)
        temp.Read(Kratos.DENSITY)

        self.properties_control.Update(temp)

    def test_PropertiesControlMapGradient(self):
        self.properties_control.Initialize()

        control_model_part = self.model_part.GetSubModelPart("test")

        with self.assertRaises(RuntimeError):
            self.properties_control.MapGradient({Kratos.DENSITY: Kratos.ContainerExpression.NodalNonHistoricalExpression(control_model_part)})

        with self.assertRaises(RuntimeError):
            self.properties_control.MapGradient({Kratos.DENSITY: KratosOA.ContainerExpression.ElementPropertiesExpression(self.model_part)})

        with self.assertRaises(RuntimeError):
            self.properties_control.MapGradient({Kratos.THICKNESS: KratosOA.ContainerExpression.ElementPropertiesExpression(control_model_part)})

        with self.assertRaises(RuntimeError):
            self.properties_control.MapGradient({
                Kratos.DENSITY: KratosOA.ContainerExpression.ElementPropertiesExpression(control_model_part),
                Kratos.THICKNESS: KratosOA.ContainerExpression.ElementPropertiesExpression(control_model_part)})

        self.properties_control.MapGradient({Kratos.DENSITY: KratosOA.ContainerExpression.ElementPropertiesExpression(control_model_part)})

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.TESTS_OUTPUTS)  # TESTS_OUTPUTS
    kratos_unittest.main()