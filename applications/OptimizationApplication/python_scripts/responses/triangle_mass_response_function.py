from typing import Optional

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartUtilities

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, _) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"TriangleMassResponseFunction instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"TriangleMassResponseFunction instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return TriangleMassResponseFunction(parameters["name"].GetString(), model, parameters["settings"])

class TriangleMassResponseFunction(ResponseFunction):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "evaluated_model_part_names"     : [
                "PLEASE_PROVIDE_A_MODEL_PART_NAME"
            ],
            "perturbation_size": 1e-6
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model = model

        evaluated_model_part_names = parameters["evaluated_model_part_names"].GetStringArray()
        if len(evaluated_model_part_names) == 0:
            raise RuntimeError(f"No model parts were provided for TriangleMassResponseFunction. [ response name = \"{self.GetName()}\"]")

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", evaluated_model_part_names, False)
        self.model_part: Optional[Kratos.ModelPart] = None
        self.perturbation_size = parameters["perturbation_size"].GetDouble()

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KratosOA.SHAPE, Kratos.DENSITY, Kratos.THICKNESS, KratosOA.CROSS_AREA]

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()

    def Check(self) -> None:
        if self.model_part is None:
            raise RuntimeError("Please call TriangleMassResponseFunction::Initialize first.")
        element = self.model_part.Elements[1]
        element.Properties[Kratos.DENSITY] = 1.0
        KratosOA.ResponseUtils.MassResponseUtils.Check(self.model_part)

    def Finalize(self) -> None:
        pass

    def GetInfluencingModelPart(self) -> Kratos.ModelPart:
        if self.model_part is None:
            raise RuntimeError("Please call TriangleMassResponseFunction::Initialize first.")
        return self.model_part

    def CalculateValue(self) -> float:
        return KratosOA.ResponseUtils.MassResponseUtils.CalculateValue(self.model_part)

    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
        # first merge all the model parts
        merged_model_part_map = ModelPartUtilities.GetMergedMap(physical_variable_collective_expressions, False)

        # now get the intersected model parts
        intersected_model_part_map = ModelPartUtilities.GetIntersectedMap(self.model_part, merged_model_part_map, False)

        # calculate the gradients
        for physical_variable, merged_model_part in merged_model_part_map.items():
            KratosOA.ResponseUtils.MassResponseUtils.CalculateGradient(
                physical_variable,
                merged_model_part,
                intersected_model_part_map[physical_variable],
                physical_variable_collective_expressions[physical_variable].GetContainerExpressions(),
                self.perturbation_size)
            
        expression = physical_variable_collective_expressions[KratosOA.SHAPE].GetContainerExpressions()[0]
        np_array = expression.Evaluate()
        np_array[0,0] = 0.0
        np_array[0,1] = 0.0
        np_array[0,2] = 0.0
        np_array[1,0] = 0.0
        np_array[1,1] = 0.0
        np_array[1,2] = 0.0
        np_array[0,2] = 0.0

        Kratos.Expression.CArrayExpressionIO.Read(expression, np_array)

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"