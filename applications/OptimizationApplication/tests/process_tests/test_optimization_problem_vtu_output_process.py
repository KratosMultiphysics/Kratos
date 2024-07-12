import KratosMultiphysics as Kratos

# Import KratosUnittest
import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.processes.optimization_problem_vtu_output_process import OptimizationProblemVtuOutputProcess
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.buffered_dict import BufferedDict
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes, SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

class TestOptimizationProblemVtuOutputProcess(kratos_unittest.TestCase):
    class DummyResponseFunction(ResponseFunction):
        def __init__(self, response_name: str, model_part: Kratos.ModelPart) -> None:
            super().__init__(response_name)
            self.model_part = model_part
        def CalculateValue(self) -> float:
            return 0.0
        def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
            pass
        def Check(self) -> None:
            pass
        def Initialize(self) -> None:
            pass
        def Finalize(self) -> None:
            pass
        def GetAnalysisModelPart(self) -> Kratos.ModelPart:
            return None
        def GetEvaluatedModelPart(self) -> Kratos.ModelPart:
            return None
        def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
            return []

    class DummyControl(Control):
        def __init__(self, control_name: str, model_part: Kratos.ModelPart) -> None:
            super().__init__(control_name)
            self.model_part = model_part
        def Check(self) -> None:
            pass
        def Initialize(self) -> None:
            pass
        def Finalize(self) -> None:
            pass
        def GetControlField(self) -> ContainerExpressionTypes:
            return None
        def GetEmptyField(self) -> ContainerExpressionTypes:
            return None
        def GetPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
            return []
        def MapGradient(self, physical_gradient_variable_container_expression_map: 'dict[SupportedSensitivityFieldVariableTypes, ContainerExpressionTypes]') -> ContainerExpressionTypes:
            return None
        def Update(self, control_field: ContainerExpressionTypes) -> bool:
            return True

    class DummyExecutionPolicy(ExecutionPolicy):
        def __init__(self, execution_policy_name: str, model_part: Kratos.ModelPart) -> None:
            super().__init__(execution_policy_name)
            self.model_part = model_part
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
        cls.model = Kratos.Model()
        cls.model_part1 = cls.model.CreateModelPart("test_1")
        cls.model_part1.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.model_part1.CreateNewNode(2, 1.0, 0.0, 0.0)
        cls.model_part1.CreateNewNode(3, 1.0, 1.0, 0.0)
        properties = cls.model_part1.CreateNewProperties(1)
        cls.model_part1.CreateNewElement("Element2D3N", 1, [1, 2, 3], properties)

        for node in cls.model_part1.Nodes:
            node.SetValue(Kratos.VELOCITY, Kratos.Array3([node.Id + 1, node.Id + 2, node.Id + 3]))

        for element in cls.model_part1.Elements:
            element.SetValue(Kratos.ACCELERATION, Kratos.Array3([element.Id + 4, element.Id + 5, element.Id + 6]))

        cls.model_part2 = cls.model.CreateModelPart("test_2")
        cls.model_part2.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.model_part2.CreateNewNode(2, 1.0, 0.0, 0.0)
        cls.model_part2.CreateNewNode(3, 1.0, 1.0, 0.0)
        cls.model_part2.CreateNewNode(4, 0.0, 1.0, 0.0)
        properties = cls.model_part2.CreateNewProperties(2)
        cls.model_part2.CreateNewElement("Element2D4N", 1, [1, 2, 3, 4], properties)

        for node in cls.model_part2.Nodes:
            node.SetValue(Kratos.VELOCITY, Kratos.Array3([node.Id + 1, node.Id + 2, node.Id + 3]))

        for element in cls.model_part2.Elements:
            element.SetValue(Kratos.ACCELERATION, Kratos.Array3([element.Id + 1, element.Id + 2, element.Id + 3]))

        cls.optimization_problem = OptimizationProblem()
        cls.components_list = []

        cls.components_list.append(TestOptimizationProblemVtuOutputProcess.DummyResponseFunction("resp_1", cls.model_part1))
        cls.components_list.append(TestOptimizationProblemVtuOutputProcess.DummyResponseFunction("resp_2", cls.model_part2))
        cls.components_list.append(TestOptimizationProblemVtuOutputProcess.DummyResponseFunction("resp_3", cls.model_part1))
        cls.components_list.append(TestOptimizationProblemVtuOutputProcess.DummyControl("control_1", cls.model_part1))
        cls.components_list.append(TestOptimizationProblemVtuOutputProcess.DummyControl("control_2", cls.model_part2))
        cls.components_list.append(TestOptimizationProblemVtuOutputProcess.DummyControl("control_3", cls.model_part1))
        cls.components_list.append(TestOptimizationProblemVtuOutputProcess.DummyExecutionPolicy("policy_1", cls.model_part1))
        cls.components_list.append(TestOptimizationProblemVtuOutputProcess.DummyExecutionPolicy("policy_2", cls.model_part2))
        cls.components_list.append(TestOptimizationProblemVtuOutputProcess.DummyExecutionPolicy("policy_3", cls.model_part1))

        for component in cls.components_list:
            cls.optimization_problem.AddComponent(component)
            ComponentDataView(component, cls.optimization_problem).SetDataBuffer(2)

    def __AddData(self, buffered_dict: BufferedDict, is_buffered_data: bool, component):
        step_v = self.optimization_problem.GetStep() + 1

        nodal_data = Kratos.Expression.NodalExpression(component.model_part)
        Kratos.Expression.VariableExpressionIO.Read(nodal_data, Kratos.VELOCITY, False)
        buffered_dict[f"{component.GetName()}_data_nodal_vel__{is_buffered_data}"] = nodal_data * step_v * 2.3

        element_data = Kratos.Expression.ElementExpression(component.model_part)
        Kratos.Expression.VariableExpressionIO.Read(element_data, Kratos.ACCELERATION)
        buffered_dict[f"{component.GetName()}_data_element_acc_{is_buffered_data}"] = element_data * step_v * 2.3

        collective_data = KratosOA.CollectiveExpression([nodal_data, element_data])
        collective_data *= 2.0
        buffered_dict[f"{component.GetName()}_collective_element_{is_buffered_data}"] = collective_data * step_v * 2.3

    def test_OptimizationProblemVtuOutputProcess(self):
        parameters = Kratos.Parameters(
            """
            {
                "file_name"                   : "<model_part_full_name>",
                "file_format"                 : "binary",
                "write_deformed_configuration": false,
                "list_of_output_components"   : ["all"],
                "output_precision"            : 7,
                "output_interval"             : 1,
                "echo_level"                  : 0
            }
            """
        )

        process = OptimizationProblemVtuOutputProcess(parameters, self.optimization_problem)
        process.ExecuteInitialize()

        # initialize unbuffered data
        for component in self.components_list:
            self.__AddData(ComponentDataView(component, self.optimization_problem).GetUnBufferedData(), False, component)

        with kratos_unittest.WorkFolderScope(".", __file__):
            number_of_steps = 10
            for _ in range(number_of_steps):
                # initialize the buffered data
                for component in self.components_list:
                    self.__AddData(ComponentDataView(component, self.optimization_problem).GetBufferedData(), True, component)

                process.PrintOutput()
                self.optimization_problem.AdvanceStep()

            process.ExecuteFinalize()

            CompareTwoFilesCheckProcess(Kratos.Parameters("""
            {
                "reference_file_name"   : "test_1_orig.vtu",
                "output_file_name"      : "Optimization_Results/test_1.vtu",
                "remove_output_file"    : true,
                "comparison_type"       : "deterministic"
            }""")).Execute()

            CompareTwoFilesCheckProcess(Kratos.Parameters("""
            {
                "reference_file_name"   : "test_2_orig.vtu",
                "output_file_name"      : "Optimization_Results/test_2.vtu",
                "remove_output_file"    : true,
                "comparison_type"       : "deterministic"
            }""")).Execute()

if __name__ == "__main__":
    kratos_unittest.main()
