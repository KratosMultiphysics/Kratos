
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.StructuralMechanicsApplication

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.OptimizationApplication.controls.shape.shape_control_new import ShapeControlNew

class TestShapeControl(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("Structure")
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3
        cls.model_part.AddNodalSolutionStepVariable(Kratos.MESH_DISPLACEMENT)
        Kratos.ModelPartIO("linear_strain_energy_test/Structure", Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(cls.model_part)

        material_settings = Kratos.Parameters("""{"Parameters": {"materials_filename": "linear_strain_energy_test/StructuralMaterials.json"}} """)
        Kratos.ReadMaterialsUtility(material_settings, cls.model)

        parameters = Kratos.Parameters("""{
            "model_part_names"                 : ["Structure.structure"],
            "mesh_motion_execution_policy_type": "IndependentAnalysisExecutionPolicy",
            "mesh_motion_module"               : "KratosMultiphysics",
            "mesh_motion_analysis_type"        : "MultistageAnalysis",
            "mesh_motion_settings"             : {
                "stages"        : [],
                "execution_list":[]
            }
        }""")

        cls.shape_control = ShapeControlNew(cls.model, parameters, None)

    @classmethod
    def tearDownClass(cls):
        with kratos_unittest.WorkFolderScope("linear_strain_energy_test", __file__):
            DeleteFileIfExisting("Structure.time")

    def test_ShapeControl(self):
        self.shape_control.Initialize()

        # run for 3 iterations
        for i in range(1, 4, 1):
            if (i > 1):
                test = Kratos.ContainerExpression.HistoricalExpression(self.model_part.GetSubModelPart("structure"))
                test.Read(Kratos.MESH_DISPLACEMENT)
                temp = Kratos.ContainerExpression.HistoricalExpression(self.model_part.GetSubModelPart("structure"))
                temp.CopyFrom(data)
                self.assertEqual(KratosOA.ContainerExpressionUtils.NormInf(test - temp), 0.0)

            data = Kratos.ContainerExpression.HistoricalExpression(self.model_part.GetSubModelPart("structure"))
            data.SetData(Kratos.Array3([i, i, i]))
            self.shape_control.Update(KratosOA.ContainerExpression.CollectiveExpressions([data]))

        self.shape_control.Finalize()

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()