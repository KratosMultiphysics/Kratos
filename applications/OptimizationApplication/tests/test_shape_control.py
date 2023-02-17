
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.StructuralMechanicsApplication

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.controls.shape_control import ShapeControl

class TestShapeControl(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.optimization_info = OptimizationInfo()
        cls.optimization_info.SetBufferSize(1)
        cls.optimization_info["step"] = 0
        cls.model_part = cls.model.CreateModelPart("Structure")
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3
        cls.model_part.AddNodalSolutionStepVariable(Kratos.MESH_DISPLACEMENT)
        Kratos.ModelPartIO("linear_element_test/Structure", Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(cls.model_part)

        material_settings = Kratos.Parameters("""{"Parameters": {"materials_filename": "linear_element_test/StructuralMaterials.json"}} """)
        Kratos.ReadMaterialsUtility(material_settings, cls.model)

        execution_policy_decorator_Executeparameters = Kratos.Parameters("""{
            "execution_policy_name"    : "test",
            "execution_policy_type"    : "IndependentAnalysisExecutionPolicy",
            "execution_policy_settings": {
                "analysis_module"  : "KratosMultiphysics",
                "analysis_type"    : "MultistageAnalysis",
                "analysis_settings": {
                    "stages": [],
                    "execution_list":[]
                }
            }
        }""")

        cls.execution_policy_decorator = ExecutionPolicyDecorator(cls.model, execution_policy_decorator_Executeparameters)
        cls.optimization_info.AddOptimizationProcess(ExecutionPolicyDecorator, cls.execution_policy_decorator.GetExecutionPolicyName(), cls.execution_policy_decorator)

        parameters = Kratos.Parameters("""{
            "model_part_names"         : ["Structure.structure"],
            "mesh_moving_analysis_name": "test"
        }""")

        cls.shape_control = ShapeControl(cls.model, parameters, cls.optimization_info)

    @classmethod
    def tearDownClass(cls):
        with kratos_unittest.WorkFolderScope("linear_element_test", __file__):
            DeleteFileIfExisting("Structure.time")

    def test_ShapeControl(self):
        self.execution_policy_decorator.ExecuteInitialize()
        self.shape_control.ExecuteInitialize()

        # run for 3 iterations
        for i in range(1, 4, 1):
            self.optimization_info.AdvanceSolutionStep()
            self.optimization_info["step"] = i
            self.execution_policy_decorator.ExecuteInitializeSolutionStep()
            self.shape_control.ExecuteInitializeSolutionStep()

            if (self.optimization_info["step"] > 1):
                test = KratosOA.HistoricalContainerVariableDataHolder(self.model_part.GetSubModelPart("structure"))
                test.ReadDataFromContainerVariable(Kratos.MESH_DISPLACEMENT)
                temp = KratosOA.HistoricalContainerVariableDataHolder(self.model_part.GetSubModelPart("structure"))
                temp.CopyDataFrom(data)
                self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.NormInf(test - temp), 0.0)

            data = KratosOA.HistoricalContainerVariableDataHolder(self.model_part.GetSubModelPart("structure"))
            data.SetDataForContainerVariable(Kratos.MESH_DISPLACEMENT, Kratos.Array3([i, i, i]))
            self.shape_control.UpdateControl(KratosOA.CollectiveVariableDataHolder([data]))

            self.execution_policy_decorator.ExecuteFinalizeSolutionStep()
            self.shape_control.ExecuteFinalizeSolutionStep()

        self.execution_policy_decorator.ExecuteFinalize()
        self.shape_control.ExecuteFinalize()

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()