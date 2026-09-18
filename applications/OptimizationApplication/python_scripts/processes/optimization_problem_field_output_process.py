from typing import Union
import abc

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import GetAllComponentFullNamesWithData
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import GetComponentHavingDataByFullName
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import TensorAdaptorData
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import CombinedTensorAdaptorData

class TensorAdaptorOutput(abc.ABC):
    @abc.abstractmethod
    def AddTensorAdaptorData(self, tensor_adaptor_data: TensorAdaptorData) -> bool:
        pass

    @abc.abstractmethod
    def WriteOutput(self):
        pass

class OptimizationProblemFieldOutputProcess(Kratos.OutputProcess):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> None:
        super().__init__()
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.model = model
        self.optimization_problem = optimization_problem
        self.parameters = parameters

        self.list_of_component_names = parameters["list_of_output_components"].GetStringArray()
        self.echo_level = parameters["echo_level"].GetInt()
        self.output_interval = parameters["output_interval"].GetInt()

        self.list_of_tensor_adaptor_outputs: 'list[TensorAdaptorOutput]' = []
        self.initialized_vtu_outputs = False
        self.last_step_written = -1

    def IsOutputStep(self) -> bool:
        return self.optimization_problem.GetStep() % self.output_interval == 0

    def PrintOutput(self) -> None:
        if not self.initialized_vtu_outputs:
            self.InitializeVtuOutputIO()
            self.initialized_vtu_outputs = True

        for tensor_adaptor_vtu_output in self.list_of_tensor_adaptor_outputs:
            tensor_adaptor_vtu_output.WriteOutput()

        self.last_step_written = self.optimization_problem.GetStep()

    def ExecuteFinalize(self) -> None:
        if self.last_step_written != self.optimization_problem.GetStep():
            # force the last step be written if it is not already written by normal interval based writing.
            self.PrintOutput()

    def InitializeVtuOutputIO(self) -> None:
        # get all the component names at the first writing point
        if len(self.list_of_component_names) == 1 and self.list_of_component_names[0] == "all":
            self.list_of_component_names = GetAllComponentFullNamesWithData(self.optimization_problem)

        list_of_components: 'list[Union[str, ResponseFunction, Control, ExecutionPolicy]]' = []
        for component_name in self.list_of_component_names:
            list_of_components.append(GetComponentHavingDataByFullName(component_name, self.optimization_problem))

        global_values_map = self.optimization_problem.GetProblemDataContainer().GetMap()
        for global_k, global_v in global_values_map.items():
             # first check whether this is part of requested list of components
            found_valid_component = False
            for component in list_of_components:
                component_data = ComponentDataView(component, self.optimization_problem)
                if global_k.startswith(component_data.GetDataPath()):
                     found_valid_component = True
                     break

            # if a valid component is found, add the tensor adaptor
            if found_valid_component:
                if isinstance(global_v, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor):
                    for i, _ in enumerate(global_v.GetTensorAdaptors()):
                        self.__AddTensorAdaptor(CombinedTensorAdaptorData(global_k, global_v, i))
                elif isinstance(global_v, Kratos.TensorAdaptors.DoubleTensorAdaptor) or  \
                     isinstance(global_v, Kratos.TensorAdaptors.IntTensorAdaptor) or \
                     isinstance(global_v, Kratos.TensorAdaptors.BoolTensorAdaptor):
                    self.__AddTensorAdaptor(TensorAdaptorData(global_k, global_v))

    def _GetModelPart(self, container) -> Kratos.ModelPart:
        def get_model_part(container, model_part: Kratos.ModelPart):
            if container in [model_part.Nodes, model_part.Conditions, model_part.Elements]:
                return model_part

            for sub_model_part_name in model_part.GetSubModelPartNames():
                root_model_part = get_model_part(container, model_part.GetSubModelPart(sub_model_part_name))
                if root_model_part is not None:
                    return root_model_part

            return None

        for model_part_name in self.model.GetModelPartNames():
            root_model_part = get_model_part(container, self.model[model_part_name])
            if root_model_part is not None:
                return root_model_part

        raise RuntimeError(f"No model part contains the provided container.")

    def __AddTensorAdaptor(self, tensor_adaptor_data: TensorAdaptorData) -> bool:
        found_vtu_output = False
        for tensor_adaptor_vtu_output in self.list_of_tensor_adaptor_outputs:
            if tensor_adaptor_vtu_output.AddTensorAdaptorData(tensor_adaptor_data):
                found_vtu_output = True
                break

        if not found_vtu_output:
            tensor_adaptor_output = self._CreateTensorAdaptorOutput(tensor_adaptor_data)
            tensor_adaptor_output.AddTensorAdaptorData(tensor_adaptor_data)
            self.list_of_tensor_adaptor_outputs.append(tensor_adaptor_output)
            if self.echo_level > 0:
                Kratos.Logger.PrintInfo(self.__class__.__name__, f"Created tensor adaptor output {tensor_adaptor_output}.")

    def _CreateTensorAdaptorOutput(self, _: TensorAdaptorData) -> TensorAdaptorOutput:
        raise NotImplementedError("_CreateTensorAdaptorOutput needs to be implemented in the derived class")

