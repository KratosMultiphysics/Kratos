from typing import Any, Union

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import GetAllComponentFullNamesWithData
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import GetComponentHavingDataByFullName
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes

def Factory(_: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Any:
    if not parameters.Has("settings"):
        raise RuntimeError(f"OptimizationProblemVtuOutputProcess instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return OptimizationProblemVtuOutputProcess(parameters["settings"], optimization_problem)

class ExpressionVtuOutput:
    def __init__(self, output_file_name_prefix: str, model_part: Kratos.ModelPart, is_initial_configuration: bool, writer_format: Kratos.VtuOutput.WriterFormat, precision: int, optimization_problem: OptimizationProblem):
        self.model_part = model_part
        self.optimization_problem = optimization_problem
        self.output_file_name_prefix = output_file_name_prefix

        self.vtu_output: Kratos.VtuOutput = Kratos.VtuOutput(model_part, is_initial_configuration, writer_format, precision)
        self.list_of_container_expression_paths: 'list[str]' = []
        self.list_of_collective_expression_paths: 'list[tuple[str, int]]' = []
        self.list_of_component_names: 'list[str]' = []

    def AddContainerExpressionPath(self, component_name: str, data_path: str):
        self._CheckIfExists(data_path)
        self.list_of_container_expression_paths.append(data_path)
        if component_name not in self.list_of_component_names:
            self.list_of_component_names.append(component_name)

    def AddCollectiveExpressionPath(self, component_name: str, data_path: str, index: int):
        self._CheckIfExists(data_path)
        self.list_of_collective_expression_paths.append([data_path, index])
        if component_name not in self.list_of_component_names:
            self.list_of_component_names.append(component_name)

    def GetModelPart(self) -> Kratos.ModelPart:
        return self.model_part

    def WriteOutput(self):
        if len(self.list_of_component_names) > 0:
            data_container = self.optimization_problem.GetProblemDataContainer()

            # now add back the new container expressions
            for container_expression_path in self.list_of_container_expression_paths:
                container_expression = data_container[container_expression_path]
                if not isinstance(container_expression, ContainerExpressionTypes):
                    raise RuntimeError(f"No container expression exists at \"{container_expression_path}\". The data container is not consistent between steps [ current data = {container_expression} ].")

                self.vtu_output.AddContainerExpression(self._GetContainerExpressionName(container_expression_path), container_expression)

            # now add the collective expressions
            for collective_expression_path, index in self.list_of_collective_expression_paths:
                collective_expression = data_container[collective_expression_path]
                if not isinstance(collective_expression, KratosOA.ContainerExpression.CollectiveExpressions):
                    raise RuntimeError(f"No collective expression exists at \"{collective_expression_path}\". The data container is not consistent between steps [ current data = {collective_expression} ].")

                self.vtu_output.AddContainerExpression(self._GetContainerExpressionName(container_expression_path), collective_expression.GetContainerExpressions()[index])

            self.vtu_output.PrintOutput(self.output_file_name_prefix + "_".join(self.list_of_component_names))

    def _CheckIfExists(self, data_path: str) -> None:
        container_expression_name = ExpressionVtuOutput._GetContainerExpressionName(data_path)

        if container_expression_name in [ExpressionVtuOutput._GetContainerExpressionName(c_path) for c_path in self.list_of_container_expression_paths]:
            raise RuntimeError(f"Failed to add the container expression. There exists already a container expression with name \"{container_expression_name}\" [ expression path = \"{data_path}\" ].")

        if container_expression_name in [ExpressionVtuOutput._GetContainerExpressionName(c_path) for c_path in self.list_of_collective_expression_paths]:
            raise RuntimeError(f"Failed to add the container expression. There exists already a collective expression with name \"{container_expression_name}\" [ expression path = \"{data_path}\" ].")

    @staticmethod
    def _GetContainerExpressionName(container_expression_path: str) -> str:
        return container_expression_path[container_expression_path.rfind("/")+1:]

class OptimizationProblemVtuOutputProcess(Kratos.OutputProcess):
    def GetDefaultParameters(self):
        return Kratos.Parameters(
            """
            {
                "custom_name_prefix"          : "SPECIFY_OUTPUT_FILE_NAME",
                "file_format"                 : "binary",
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
        self.output_name_prefix = parameters["custom_name_prefix"].GetString()
        self.write_deformed_configuration = parameters["write_deformed_configuration"].GetBool()
        self.output_precision = parameters["output_precision"].GetInt()
        file_format = parameters["file_format"].GetString()
        if file_format == "ascii":
            self.writer_format = Kratos.VtuOutput.ASCII
        elif file_format == "binary":
            self.writer_format = Kratos.VtuOutput.BINARY
        else:
            raise RuntimeError(f"Only supports \"ascii\" and \"binary\" file_format. [ provided file_format = \"{file_format}\" ].")

        self.list_of_components: 'list[Union[str, ResponseFunction, Control, ExecutionPolicy]]' = []
        list_of_component_names = parameters["list_of_output_components"].GetStringArray()
        if len(list_of_component_names) == 1 and list_of_component_names[0] == "all":
            list_of_component_names = GetAllComponentFullNamesWithData(optimization_problem)

        for component_name in list_of_component_names:
            self.list_of_components.append(GetComponentHavingDataByFullName(component_name, optimization_problem))

        self.list_of_expresson_vtu_outputs: 'list[ExpressionVtuOutput]' = []
        self.initialized_vtu_outputs = False

    def PrintOutput(self) -> None:
        if not self.initialized_vtu_outputs:
            self.InitializeVtuOutputIO()
            self.initialized_vtu_outputs = True

        for expression_vtu_output in self.list_of_expresson_vtu_outputs:
            expression_vtu_output.WriteOutput()

    def InitializeVtuOutputIO(self) -> None:
        global_values_map = self.optimization_problem.GetProblemDataContainer().GetMap()
        for global_k, global_v in global_values_map.items():
             # first check whether this is part of requested list of components
            found_valid_component = False
            for component in self.list_of_components:
                component_data = ComponentDataView(component, self.optimization_problem)
                if global_k.startswith(component_data.GetDataPath()):
                     found_valid_component = True
                     break

            # if a valid component is found, add the expression
            if found_valid_component:
                if isinstance(global_v, ContainerExpressionTypes):
                    model_part = global_v.GetModelPart()
                    self._AddContainerExpression(component_data, global_k, model_part)
                elif isinstance(global_v, KratosOA.ContainerExpression.CollectiveExpressions):
                    for i, container_expression in enumerate(global_v.GetContainerExpressions()):
                        model_part = container_expression.GetModelPart()
                        self._AddContainerExpression(component_data, global_k, model_part, i)

    def _AddContainerExpression(self, component_data: ComponentDataView, data_path: str, model_part: Kratos.ModelPart, index: int = -1):
        found_vtu_output = False
        for expression_vtu_output in self.list_of_expresson_vtu_outputs:
            if model_part == expression_vtu_output.GetModelPart():
                found_vtu_output = True
                if index == -1:
                    expression_vtu_output.AddContainerExpressionPath(component_data.GetComponentName(), data_path)
                else:
                    expression_vtu_output.AddCollectiveExpressionPath(component_data.GetComponentName(), data_path, index)
                break

        if not found_vtu_output:
            expression_vtu_output = ExpressionVtuOutput(self.output_name_prefix, model_part, not self.write_deformed_configuration, self.writer_format, self.output_precision, self.optimization_problem)
            if index == -1:
                expression_vtu_output.AddContainerExpressionPath(component_data.GetComponentName(), data_path)
            else:
                expression_vtu_output.AddCollectiveExpressionPath(component_data.GetComponentName(), data_path, index)

            self.list_of_expresson_vtu_outputs.append(expression_vtu_output)
            if self.echo_level > 0:
                Kratos.Logger.PrintInfo(self.__class__.__name__, f"Created expression vtu output for {model_part.FullName()}.")