import KratosMultiphysics as Kratos

# Import KratosUnittest
import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.processes.optimization_problem_vtu_output_process import OptimizationProblemVtuOutputProcess
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.buffered_dict import BufferedDict
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess


class TestOptimizationProblemVtuOutputProcess(kratos_unittest.TestCase):
    class DummyResponseFunction(ResponseFunction):
        def __init__(self, response_name: str, model_part: Kratos.ModelPart) -> None:
            super().__init__(response_name)
            self.model_part = model_part
        def CalculateValue(self) -> float:
            return 0.0
        def CalculateGradient(self, physical_variable_combined_tensor_adaptors: 'dict[SupportedSensitivityFieldVariableTypes, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]') -> None:
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
        def __init__(self, control_name: str, model_part: Kratos.ModelPart) -> None:
            super().__init__(control_name)
            self.model_part = model_part
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
        def MapGradient(self, physical_gradient_variable_tensor_adaptor_map: 'dict[SupportedSensitivityFieldVariableTypes, Kratos.TensorAdaptors.DoubleTensorAdaptor]') -> Kratos.TensorAdaptors.DoubleTensorAdaptor:
            return None
        def Update(self, control_field: Kratos.TensorAdaptors.DoubleTensorAdaptor) -> bool:
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
        # purposefully ignored to always overwrite to the same step = 0 file.
        # self.model_part1.ProcessInfo[Kratos.STEP] = step_v
        # self.model_part2.ProcessInfo[Kratos.STEP] = step_v
        self.model_part1.ProcessInfo[Kratos.TIME] = step_v
        self.model_part2.ProcessInfo[Kratos.TIME] = step_v

        nodal_data = Kratos.TensorAdaptors.VariableTensorAdaptor(component.model_part.Nodes, Kratos.VELOCITY)
        nodal_data.CollectData()
        nodal_data.data[:] *= step_v * 2.3
        buffered_dict[f"{component.GetName()}_data_nodal_vel__{is_buffered_data}"] = nodal_data

        element_data = Kratos.TensorAdaptors.VariableTensorAdaptor(component.model_part.Elements, Kratos.ACCELERATION)
        element_data.CollectData()
        element_data.data[:] *= step_v * 2.3
        buffered_dict[f"{component.GetName()}_data_element_acc_{is_buffered_data}"] = element_data

        combined_data = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor([nodal_data, element_data])
        combined_data.CollectData()
        combined_data.data[:] *= 2.0 * step_v * 2.3
        Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(combined_data, perform_store_data_recursively=False, copy=False).StoreData()
        buffered_dict[f"{component.GetName()}_combined_element_{is_buffered_data}"] = combined_data

    def test_OptimizationProblemVtuOutputProcess(self):
        parameters = Kratos.Parameters(
            """
            {
                "file_name"                   : "<model_part_full_name>",
                "file_format"                 : "ascii",
                "write_deformed_configuration": false,
                "list_of_output_components"   : ["all"],
                "output_precision"            : 7,
                "output_interval"             : 1,
                "echo_level"                  : 0
            }
            """
        )

        process = OptimizationProblemVtuOutputProcess(self.model, parameters, self.optimization_problem)
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
                "output_file_name"      : "Optimization_Results/test_1/test_1_elements_0.vtu",
                "remove_output_file"    : true,
                "comparison_type"       : "deterministic"
            }""")).Execute()

            CompareTwoFilesCheckProcess(Kratos.Parameters("""
            {
                "reference_file_name"   : "test_2_orig.vtu",
                "output_file_name"      : "Optimization_Results/test_2/test_2_elements_0.vtu",
                "remove_output_file"    : true,
                "comparison_type"       : "deterministic"
            }""")).Execute()

    @classmethod
    def tearDownClass(cls):
        with kratos_unittest.WorkFolderScope(".", __file__):
            kratos_utilities.DeleteDirectoryIfExisting("Optimization_Results")

if __name__ == "__main__":
    kratos_unittest.main()
