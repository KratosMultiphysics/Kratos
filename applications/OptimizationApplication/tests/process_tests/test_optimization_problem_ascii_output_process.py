
import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.processes.optimization_problem_ascii_output_process import OptimizationProblemAsciiOutputProcess
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.buffered_dict import BufferedDict
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

class TestOptimizationProblemAsciiOutputProcess(kratos_unittest.TestCase):
    class DummyResponseFunction(ResponseFunction):
        def __init__(self, response_name: str) -> None:
            super().__init__(response_name)
        def CalculateValue(self) -> float:
            return 0.0
        def CalculateGradient(self, physical_variable_combined_tensor_adaptor: 'dict[SupportedSensitivityFieldVariableTypes, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]') -> None:
            pass
        def Check(self) -> None:
            pass
        def Initialize(self) -> None:
            pass
        def Finalize(self) -> None:
            pass
        def GetInfluencingModelPart(self) -> Kratos.ModelPart:
            return None
        def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
            return []

    class DummyControl(Control):
        def __init__(self, control_name: str) -> None:
            super().__init__(control_name)
        def Check(self) -> None:
            pass
        def Initialize(self) -> None:
            pass
        def Finalize(self) -> None:
            pass
        def GetControlField(self) -> Kratos.TensorAdaptors.DoubleTensorAdaptor:
            return None
        def GetEmptyField(self) -> Kratos.TensorAdaptors.DoubleTensorAdaptor:
            return None
        def GetPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
            return []
        def MapGradient(self, physical_gradient_variable_container_expression_map: 'dict[SupportedSensitivityFieldVariableTypes, Kratos.TensorAdaptors.DoubleTensorAdaptor]') -> Kratos.TensorAdaptors.DoubleTensorAdaptor:
            return None
        def Update(self, control_field: Kratos.TensorAdaptors.DoubleTensorAdaptor) -> bool:
            return True

    class DummyExecutionPolicy(ExecutionPolicy):
        def __init__(self, execution_policy_name: str) -> None:
            super().__init__(execution_policy_name)
        def Check(self) -> None:
            pass
        def Initialize(self) -> None:
            pass
        def Finalize(self) -> None:
            pass
        def GetAnalysisModelPart(self) -> Kratos.ModelPart:
            return None
        def Execute(self) -> None:
            pass

    @classmethod
    def setUpClass(cls) -> None:
        cls.optimization_problem = OptimizationProblem()
        cls.components_list = []

        cls.components_list.append(TestOptimizationProblemAsciiOutputProcess.DummyResponseFunction("resp_1"))
        cls.components_list.append(TestOptimizationProblemAsciiOutputProcess.DummyResponseFunction("resp_2"))
        cls.components_list.append(TestOptimizationProblemAsciiOutputProcess.DummyResponseFunction("resp_3"))
        cls.components_list.append(TestOptimizationProblemAsciiOutputProcess.DummyControl("control_1"))
        cls.components_list.append(TestOptimizationProblemAsciiOutputProcess.DummyControl("control_2"))
        cls.components_list.append(TestOptimizationProblemAsciiOutputProcess.DummyControl("control_3"))
        cls.components_list.append(TestOptimizationProblemAsciiOutputProcess.DummyExecutionPolicy("policy_1"))
        cls.components_list.append(TestOptimizationProblemAsciiOutputProcess.DummyExecutionPolicy("policy_2"))
        cls.components_list.append(TestOptimizationProblemAsciiOutputProcess.DummyExecutionPolicy("policy_3"))

        for component in cls.components_list:
            cls.optimization_problem.AddComponent(component)
            ComponentDataView(component, cls.optimization_problem).SetDataBuffer(2)

        cls.components_list.append("algorithm")
        ComponentDataView("algorithm", cls.optimization_problem).SetDataBuffer(2)


    def __AddData(self, buffered_dict: BufferedDict):
        step_v = self.optimization_problem.GetStep() + 1
        buffered_dict[f"v_float"] = 2.3 * step_v
        buffered_dict[f"v_int"] = 2 * step_v
        buffered_dict[f"v_bool"] = bool(step_v % 2)
        buffered_dict[f"v_str"] = f"lengthy_str_{step_v}"

    def test_OptimizationProblemAsciiOutputProcess(self):
        parameters = Kratos.Parameters(
            """
            {
                "output_file_name"         : "dummy_output_test.csv",
                "write_kratos_version"     : false,
                "write_time_stamp"         : false,
                "list_of_output_components": ["all"],
                "format_info": {
                    "int_length"     : 7,
                    "float_precision": 9,
                    "bool_values"    : ["no", "yes"],
                    "string_length"  : 14
                }
            }
            """
        )

        process = OptimizationProblemAsciiOutputProcess(parameters, self.optimization_problem)

        # initialize unbuffered data
        for component in self.components_list:
            self.__AddData(ComponentDataView(component, self.optimization_problem).GetUnBufferedData())

        with kratos_unittest.WorkFolderScope(".", __file__):
            number_of_steps = 10
            for _ in range(number_of_steps):
                # initialize the buffered data
                for component in self.components_list:
                    self.__AddData(ComponentDataView(component, self.optimization_problem).GetBufferedData())

                process.PrintOutput()
                self.optimization_problem.AdvanceStep()

            process.ExecuteFinalize()

            CompareTwoFilesCheckProcess(Kratos.Parameters("""
            {
                "reference_file_name"   : "dummy_output_test_orig.csv",
                "output_file_name"      : "dummy_output_test.csv",
                "remove_output_file"    : true,
                "comparison_type"       : "deterministic"
            }""")).Execute()

if __name__ == "__main__":
    kratos_unittest.main()