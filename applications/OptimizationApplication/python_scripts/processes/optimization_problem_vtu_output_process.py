from typing import Any, Union
from abc import ABC, abstractmethod
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.kratos_utilities as kratos_utils
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import GetAllComponentFullNamesWithData
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import GetComponentHavingDataByFullName
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes

def Factory(_: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Any:
    if not parameters.Has("settings"):
        raise RuntimeError(f"OptimizationProblemVtuOutputProcess instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return OptimizationProblemVtuOutputProcess(parameters["settings"], optimization_problem)

class ExpressionData(ABC):
    def __init__(self, expression_path: str, container_expression: ContainerExpressionTypes) -> None:
        self.expression_path = expression_path
        self.container_expression_type = type(container_expression)
        self.model_part = container_expression.GetModelPart()

    def GetContainerExpressionType(self) -> Any:
        return self.container_expression_type

    def GetContainerExpressionName(self) -> str:
        return self.expression_path[self.expression_path.rfind("/")+1:]

    def GetExpressionPath(self) -> str:
        return self.expression_path

    def GetModelPart(self) -> Kratos.ModelPart:
        return self.model_part

    def _IsValidContainerExpression(self, container_expression: ContainerExpressionTypes) -> bool:
        return isinstance(container_expression, self.container_expression_type) and self.model_part == container_expression.GetModelPart()

    @abstractmethod
    def GetContainerExpression(self, optiization_problem: OptimizationProblem) -> ContainerExpressionTypes:
        pass

class ContainerExpressionData(ExpressionData):
    def __inif__(self, expression_path: str, container_expression: ContainerExpressionTypes):
        super().__init__(expression_path, container_expression)

    def GetContainerExpression(self, optiization_problem: OptimizationProblem) -> ContainerExpressionTypes:
        data = optiization_problem.GetProblemDataContainer()[self.expression_path]
        if not self._IsValidContainerExpression(data):
            raise RuntimeError(f"The data at \"{self.expression_path}\" is not valid. The original data found at this location is of type {self.container_expression_type.__name__} with {self.model_part.FullName()}. Found data = {data}")
        return data

class CollectiveExpressionData(ExpressionData):
    def __init__(self, collective_expression_path: str, collective_expression: KratosOA.CollectiveExpression, container_expression_index: int):
        super().__init__(collective_expression_path, collective_expression.GetContainerExpressions()[container_expression_index])
        self.container_expression_index = container_expression_index

    def GetContainerExpression(self, optiization_problem: OptimizationProblem) -> ContainerExpressionTypes:
        data = optiization_problem.GetProblemDataContainer()[self.expression_path]
        if not isinstance(data, KratosOA.CollectiveExpression):
            raise RuntimeError(f"The data type at \"{self.expression_path}\" changed from {KratosOA.CollectiveExpression.__name__} to {type(data).__class__.__name__}.")
        data = data.GetContainerExpressions()[self.container_expression_index]
        if not self._IsValidContainerExpression(data):
            raise RuntimeError(f"The data at \"{self.expression_path}\" is not valid. The original data found at this location is of type {self.container_expression_type.__name__} with {self.model_part.FullName()}. Found data = {data}")
        return data

class ExpressionVtuOutput:
    def __init__(self, parameters: 'dict[str, Any]', model_part: Kratos.ModelPart, optimization_problem: OptimizationProblem):
        self.model_part = model_part
        self.optimization_problem = optimization_problem
        self.output_file_name_prefix = parameters["output_file_name_prefix"]

        if parameters["save_output_files_in_folder"]:
            self.output_path = Path(parameters["output_path"])
            if not self.model_part.ProcessInfo[Kratos.IS_RESTARTED]:
                kratos_utils.DeleteDirectoryIfExisting(str(self.output_path))
            self.model_part.GetCommunicator().GetDataCommunicator().Barrier()
            # now create the output path
            Kratos.FilesystemExtensions.MPISafeCreateDirectories(str(self.output_path))
        else:
            self.output_path = Path(".")

        self.vtu_output: Kratos.VtuOutput = Kratos.VtuOutput(model_part, parameters["is_initial_configuration"], parameters["writer_format"], parameters["precision"])
        self.dict_of_expression_data: 'dict[Any, dict[str, ExpressionData]]' = {}

        # vtu output gives priority to elements over conditions if both are present.
        # nodal container expression is allowed in any case. hence that is added first
        self.dict_of_expression_data[Kratos.Expression.NodalExpression]  = {}

        # now check for elements
        communicator: Kratos.Communicator = self.model_part.GetCommunicator()
        if communicator.GlobalNumberOfElements() > 0:
            self.dict_of_expression_data[Kratos.Expression.ElementExpression] = {}
        elif communicator.GlobalNumberOfConditions() > 0:
            self.dict_of_expression_data[Kratos.Expression.ConditionExpression] = {}

    def AddExpressionData(self, expression_data: ExpressionData) -> bool:
        if expression_data.GetModelPart() == self.vtu_output.GetModelPart():
            current_expression_type = expression_data.GetContainerExpressionType()
            if current_expression_type not in self.dict_of_expression_data.keys():
                raise RuntimeError(f"The {current_expression_type.__name__} is not supported to be written in vtu format for {self.model_part.FullName()}")

            dict_of_current_expression_type = self.dict_of_expression_data[current_expression_type]

            current_expression_name = expression_data.GetContainerExpressionName()
            if current_expression_name not in dict_of_current_expression_type.keys():
                dict_of_current_expression_type[current_expression_name] = expression_data
            else:
                raise RuntimeError(f"Trying to add a duplicate expression with the name = \"{current_expression_name}\" [ The original data path = \"{dict_of_current_expression_type[current_expression_name].GetExpressionPath()}\", current data path = \"{expression_data.GetExpressionPath()}\"].")

            return True

        return False

    def WriteOutput(self):
        for current_expression_name_data_pair in self.dict_of_expression_data.values():
            for expression_data in current_expression_name_data_pair.values():
                self.vtu_output.AddContainerExpression(expression_data.GetContainerExpressionName(), expression_data.GetContainerExpression(self.optimization_problem))

        output_file_name = self.output_file_name_prefix
        output_file_name = output_file_name.replace("<model_part_full_name>", self.model_part.FullName())
        output_file_name = output_file_name.replace("<model_part_name>", self.model_part.Name)
        output_file_name = output_file_name.replace("<step>", str(self.optimization_problem.GetStep()))
        self.vtu_output.PrintOutput(str(self.output_path /output_file_name))

class OptimizationProblemVtuOutputProcess(Kratos.OutputProcess):
    def GetDefaultParameters(self):
        return Kratos.Parameters(
            """
            {
                "file_name"                   : "<model_part_full_name>_<step>",
                "file_format"                 : "binary",
                "output_path"                 : "Optimization_Results",
                "save_output_files_in_folder" : true,
                "write_deformed_configuration": false,
                "list_of_output_components"   : ["all"],
                "output_precision"            : 7,
                "output_interval"             : 1,
                "echo_level"                  : 0
            }
            """
        )

    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> None:
        super().__init__()

        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.optimization_problem = optimization_problem

        self.echo_level = parameters["echo_level"].GetInt()
        self.file_name = parameters["file_name"].GetString()
        self.output_path = parameters["output_path"].GetString()
        self.save_output_files_in_folder = parameters["save_output_files_in_folder"].GetBool()
        self.write_deformed_configuration = parameters["write_deformed_configuration"].GetBool()
        self.output_precision = parameters["output_precision"].GetInt()
        file_format = parameters["file_format"].GetString()
        if file_format == "ascii":
            self.writer_format = Kratos.VtuOutput.ASCII
        elif file_format == "binary":
            self.writer_format = Kratos.VtuOutput.BINARY
        else:
            raise RuntimeError(f"Only supports \"ascii\" and \"binary\" file_format. [ provided file_format = \"{file_format}\" ].")

        self.list_of_component_names = parameters["list_of_output_components"].GetStringArray()
        self.list_of_expresson_vtu_outputs: 'list[ExpressionVtuOutput]' = []
        self.initialized_vtu_outputs = False

    def PrintOutput(self) -> None:
        if not self.initialized_vtu_outputs:
            self.InitializeVtuOutputIO()
            self.initialized_vtu_outputs = True

        for expression_vtu_output in self.list_of_expresson_vtu_outputs:
            expression_vtu_output.WriteOutput()

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

            # if a valid component is found, add the expression
            if found_valid_component:
                if isinstance(global_v, Kratos.Expression.NodalExpression) or \
                   isinstance(global_v, Kratos.Expression.ConditionExpression) or \
                   isinstance(global_v, Kratos.Expression.ElementExpression):
                    self._AddContainerExpression(ContainerExpressionData(global_k, global_v))
                elif isinstance(global_v, KratosOA.CollectiveExpression):
                    for i, _ in enumerate(global_v.GetContainerExpressions()):
                        self._AddContainerExpression(CollectiveExpressionData(global_k, global_v, i))

    def _AddContainerExpression(self, expression_data: ExpressionData):
        found_vtu_output = False
        for expression_vtu_output in self.list_of_expresson_vtu_outputs:
            if expression_vtu_output.AddExpressionData(expression_data):
                found_vtu_output = True
                break

        if not found_vtu_output:
            vtu_parameters = {
                "output_file_name_prefix": self.file_name,
                "is_initial_configuration": not self.write_deformed_configuration,
                "writer_format": self.writer_format,
                "precision": self.output_precision,
                "save_output_files_in_folder": self.save_output_files_in_folder,
                "output_path": self.output_path
            }

            expression_vtu_output = ExpressionVtuOutput(vtu_parameters, expression_data.GetModelPart(), self.optimization_problem)
            expression_vtu_output.AddExpressionData(expression_data)
            self.list_of_expresson_vtu_outputs.append(expression_vtu_output)
            if self.echo_level > 0:
                Kratos.Logger.PrintInfo(self.__class__.__name__, f"Created expression vtu output for {expression_data.GetModelPart().FullName()}.")