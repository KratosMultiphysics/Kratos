from typing import Optional

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartUtilities

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, _) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"GeometricCentroidDeviationResponseFunction instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"GeometricCentroidDeviationResponseFunction instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return GeometricCentroidDeviationResponseFunction(parameters["name"].GetString(), model, parameters["settings"])

class GeometricCentroidDeviationResponseFunction(ResponseFunction):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "evaluated_model_part_names"     : [
                "PLEASE_PROVIDE_A_MODEL_PART_NAME"
            ]
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model = model

        evaluated_model_part_names = parameters["evaluated_model_part_names"].GetStringArray()
        if len(evaluated_model_part_names) == 0:
            raise RuntimeError(f"No model parts were provided for GeometricCentroidDeviationResponseFunction. [ response name = \"{self.GetName()}\"]")

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", evaluated_model_part_names, False)
        self.model_part: Optional[Kratos.ModelPart] = None

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KratosOA.SHAPE]

    def GetEvaluatedModelPart(self) -> Kratos.ModelPart:
        if self.model_part is None:
            raise RuntimeError("Please call GeometricCentroidDeviationResponseFunction::Initialize first.")
        return self.model_part

    def GetAnalysisModelPart(self) -> None:
        return None

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()

        self.model_part_center = Kratos.Array3(0.0)
        number_of_nodes = self.model_part.GetCommunicator().GlobalNumberOfNodes()

        if number_of_nodes == 0:
            raise RuntimeError(f"The provided model part {self.model_part.FullName()} does not contain any nodes.")

        for node in self.model_part.Nodes:
            self.model_part_center[0] += node.X
            self.model_part_center[1] += node.Y
            self.model_part_center[2] += node.Z

        self.model_part_center = self.model_part.GetCommunicator().GetDataCommunicator().SumAll(self.model_part_center)
        self.model_part_center /= number_of_nodes

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def CalculateValue(self) -> float:
        average_location = Kratos.Array3(0.0)
        for node in self.model_part.Nodes:
            average_location[0] += node.X
            average_location[1] += node.Y
            average_location[2] += node.Z

        average_location = self.model_part.GetCommunicator().GetDataCommunicator().SumAll(average_location)

        number_of_nodes = self.model_part.GetCommunicator().GlobalNumberOfNodes()
        self.value_array = (average_location / number_of_nodes  - self.model_part_center)
        return self.value_array[0] ** 2 + self.value_array[1] ** 2 + self.value_array[2] ** 2

    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
        # first merge all the model parts
        merged_model_part_map = ModelPartUtilities.GetMergedMap(physical_variable_collective_expressions, False)

        # now get the intersected model parts
        intersected_model_part_map = ModelPartUtilities.GetIntersectedMap(self.model_part, merged_model_part_map, False)

        number_of_nodes = self.model_part.GetCommunicator().GlobalNumberOfNodes()

        # calculate the gradients
        for physical_variable, merged_model_part in merged_model_part_map.items():
            if physical_variable == KratosOA.SHAPE:
                Kratos.VariableUtils().SetNonHistoricalVariableToZero(Kratos.SHAPE_SENSITIVITY, merged_model_part.Nodes)
                Kratos.VariableUtils().SetNonHistoricalVariable(Kratos.SHAPE_SENSITIVITY, 2.0 * self.value_array / number_of_nodes, intersected_model_part_map[physical_variable].Nodes)
                for container_expression in physical_variable_collective_expressions[physical_variable].GetContainerExpressions():
                    if isinstance(container_expression, Kratos.Expression.NodalExpression):
                        Kratos.Expression.VariableExpressionIO.Read(container_expression, Kratos.SHAPE_SENSITIVITY, False)
                    else:
                        raise RuntimeError(f"Requesting sensitivity w.r.t. SHAPE for a Container expression which is not a NodalExpression. [ Requested container expression = {container_expression} ].")
            else:
                raise RuntimeError(f"Unsupported sensitivity w.r.t. {physical_variable.Name()} requested. Followings are supported sensitivity variables:\n\tSHAPE")
