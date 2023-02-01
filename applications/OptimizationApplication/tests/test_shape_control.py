
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.StructuralMechanicsApplication

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_wrapper import ExecutionPolicyWrapper
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

        execution_policy_wrapper_parameters = Kratos.Parameters("""{
            "name"                     : "test",
            "execution_policy_settings": {
                "type"    : "IndependentAnalysisExecutionPolicy",
                "settings": {
                    "analysis_settings": {
                        "module"  : "KratosMultiphysics",
                        "type"    : "MultistageAnalysis",
                        "settings": {
                            "stages": [],
                            "execution_list":[]
                        }
                    }
                }
            }
        }""")

        cls.execution_policy_wrapper = ExecutionPolicyWrapper(cls.model, execution_policy_wrapper_parameters)
        cls.optimization_info.AddOptimizationRoutine(ExecutionPolicyWrapper, cls.execution_policy_wrapper.GetName(), cls.execution_policy_wrapper)

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
        self.execution_policy_wrapper.Initialize()
        self.shape_control.Initialize()

        # run for 3 iterations
        for i in range(1, 4, 1):
            self.optimization_info.AdvanceSolutionStep()
            self.optimization_info["step"] = i
            self.execution_policy_wrapper.InitializeSolutionStep()
            self.shape_control.InitializeSolutionStep()

            if (self.optimization_info["step"] > 1):
                test = KratosOA.HistoricalContainerVariableDataHolder(self.model_part.GetSubModelPart("structure"))
                test.ReadDataFromContainerVariable(Kratos.MESH_DISPLACEMENT)
                temp = KratosOA.HistoricalContainerVariableDataHolder(self.model_part.GetSubModelPart("structure"))
                temp.CopyDataFrom(data)
                self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.NormInf(test - temp), 0.0)

            data = self.shape_control.CreateContainerVariableDataHolder(self.model_part.GetSubModelPart("structure"))
            data.SetDataForContainerVariable(Kratos.MESH_DISPLACEMENT, Kratos.Array3([i, i, i]))
            self.shape_control.UpdateControl(data)

            self.execution_policy_wrapper.FinalizeSolutionStep()
            self.shape_control.FinalizeSolutionStep()

        self.execution_policy_wrapper.Finalize()
        self.shape_control.Finalize()

    def test_GetModelPart(self):
        self.assertEqual(self.model_part.GetSubModelPart("structure"), self.shape_control.GetModelParts()[0])

    def test_GetContainerType(self):
        self.assertTrue(isinstance(self.shape_control.CreateContainerVariableDataHolder(self.model_part), KratosOA.NodalContainerVariableDataHolder))

    def test_GetControlSensitivityVariable(self):
        self.assertEqual(self.shape_control.GetControlSensitivityVariable(), Kratos.SHAPE_SENSITIVITY)

    def test_GetControlUpdateVariable(self):
        self.assertEqual(self.shape_control.GetControlUpdateVariable(), KratosOA.VECTOR_CONTROL_UPDATE)

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()