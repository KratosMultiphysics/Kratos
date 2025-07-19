import numpy
from typing import Optional
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"ExternalResponseFunction instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"ExternalResponseFunction instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return ExternalResponseFunction(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)

class ExternalResponseFunction(ResponseFunction):
    def __init__(self, response_name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> None:
        super().__init__(response_name)

        default_parameters = Kratos.Parameters("""{
            "evaluated_model_part_names": []
        }""")
        parameters.ValidateAndAssignDefaults(default_parameters)

        self.model = model
        self.optimization_problem = optimization_problem

        evaluated_model_part_names = parameters["evaluated_model_part_names"].GetStringArray()
        if len(evaluated_model_part_names) == 0:
            raise RuntimeError(f"No model parts were provided for MassResponseFunction. [ response name = \"{self.GetName()}\"]")
        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", evaluated_model_part_names, False)
        self.model_part: Optional[Kratos.ModelPart] = None

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KratosOA.CUSTOM_DESIGN_VARIABLE]

    def GetInfluencingModelPart(self) -> Kratos.ModelPart:
        return self.model_part

    def _GetDesignVariableValues(self) -> numpy.ndarray:
        # get the values from the model part nodal non-historical container
        exp = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(exp, KratosOA.CUSTOM_DESIGN_VARIABLE, False)

        # convert to numpy
        return exp.Evaluate()

    def _UpdateGradientsFromNumpy(self, numpy_values: numpy.ndarray, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
        # assign the sensitivities to the physical_variable_collective_expressions
        Kratos.Expression.CArrayExpressionIO.Read(physical_variable_collective_expressions[KratosOA.CUSTOM_DESIGN_VARIABLE].GetContainerExpressions()[0], numpy_values)

    def CalculateValue(self) -> float:
        numpy_values = self._GetDesignVariableValues()

        # do the update of the model here with the new values.
        # there after compute the new response value.

        return float((numpy_values[0]-2) ** 2 + (numpy_values[1]-3) ** 2 + (numpy_values[2]-4) ** 2 + (numpy_values[3]-5))**2

    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
        numpy_values = self._GetDesignVariableValues()

        # compute the gradients
        df_dx0 = 2 * (numpy_values[0]-2)
        df_dx1 = 2 * (numpy_values[1]-3)
        df_dx2 = 2 * (numpy_values[2]-4)
        df_dx3 = 2 * (numpy_values[3]-5)

        # create the numpy array of sensitivities
        numpy_exp = numpy.array([df_dx0, df_dx1, df_dx2, df_dx3])

        self._UpdateGradientsFromNumpy(numpy_exp, physical_variable_collective_expressions)
