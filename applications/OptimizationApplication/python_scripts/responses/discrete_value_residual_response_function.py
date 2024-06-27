from typing import Optional

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, _) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"DiscreteValueResidualResponseFunction instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"DiscreteValueResidualResponseFunction instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return DiscreteValueResidualResponseFunction(parameters["name"].GetString(), model, parameters["settings"])

class DiscreteValueResidualResponseFunction(ResponseFunction):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "evaluated_model_part_names"     : [
                "PLEASE_PROVIDE_A_MODEL_PART_NAME"
            ],
            "container_type"         : "node_historical",
            "variable_name"          : "",
            "residual_type"          : "exact",
            "list_of_discrete_values": [0.0]
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model = model
        self.variable: SupportedSensitivityFieldVariableTypes = Kratos.KratosGlobals.GetVariable(parameters["variable_name"].GetString())
        if not isinstance(self.variable, Kratos.DoubleVariable):
            raise RuntimeError(f"Only supports double variables. Provided variable is \"{self.variable.Name()}\".")

        self.list_of_discrete_values = parameters["list_of_discrete_values"].GetVector()

        self.residual_type = parameters["residual_type"].GetString()
        if self.residual_type  not in ["exact", "logarithm"]:
            raise RuntimeError(f"Unsupported residual_type = \"{self.residual_type}\" requested. Followings are supported:\n\texact\n\tlogarithm")

        container_type = parameters["container_type"].GetString()
        if container_type == "node_historical" or container_type == "node_non_historical":
            self.expression_getter = lambda model_part : Kratos.Expression.NodalExpression(model_part)
            self.expression_reader = lambda exp: Kratos.Expression.VariableExpressionIO.Read(exp, self.variable, container_type == "node_historical")
        elif container_type == "condition" or container_type == "condition_properties":
            self.expression_getter = lambda model_part : Kratos.Expression.ConditionExpression(model_part)
            if container_type == "condition":
                self.expression_reader = lambda exp: Kratos.Expression.VariableExpressionIO.Read(exp, self.variable)
            else:
                self.expression_reader = lambda exp: KratosOA.PropertiesVariableExpressionIO(exp, self.variable)
        elif container_type == "element" or container_type == "element_properties":
            self.expression_getter = lambda model_part : Kratos.Expression.ElementExpression(model_part)
            if container_type == "element":
                self.expression_reader = lambda exp: Kratos.Expression.VariableExpressionIO.Read(exp, self.variable)
            else:
                self.expression_reader = lambda exp: KratosOA.PropertiesVariableExpressionIO(exp, self.variable)
        else:
            raise RuntimeError(f"Unsupported container_type = \"{container_type}\" requested. Followings are supported:\n\tnode_historical\n\tnode_non_historical\n\tcondition\n\tcondition_properties\n\telement\n\telement_properties")

        evaluated_model_part_names = parameters["evaluated_model_part_names"].GetStringArray()
        if len(evaluated_model_part_names) == 0:
            raise RuntimeError(f"No model parts were provided for DiscreteValueResidualResponseFunction. [ response name = \"{self.GetName()}\"]")

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", evaluated_model_part_names, False)
        self.model_part: Optional[Kratos.ModelPart] = None

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [self.variable]

    def GetEvaluatedModelPart(self) -> Kratos.ModelPart:
        if self.model_part is None:
            raise RuntimeError("Please call DiscreteValueResidualResponseFunction::Initialize first.")
        return self.model_part

    def GetAnalysisModelPart(self) -> None:
        return None

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def CalculateValue(self) -> float:
        exp = self.expression_getter(self.model_part)
        self.expression_reader(exp)

        resultant = exp.Clone()

        if self.residual_type == "exact":
            Kratos.Expression.LiteralExpressionIO.SetData(resultant, 1.0)
        elif self.residual_type == "logarithm":
            Kratos.Expression.LiteralExpressionIO.SetData(resultant, 0.0)

        for value in self.list_of_discrete_values:
            if self.residual_type == "exact":
                resultant *= (exp - value) ** 2
            elif self.residual_type == "logarithm":
                resultant += Kratos.Expression.Utils.Log((exp - value) ** 2)

        return Kratos.Expression.Utils.Sum(resultant)

    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
        values = self.expression_getter(self.model_part)
        self.expression_reader(values)

        # calculate the gradients
        for physical_variable, collective_expression in physical_variable_collective_expressions.items():
            if physical_variable == self.variable:
                expressions = collective_expression.GetContainerExpressions()

                # initialize the current expression
                for exp in expressions:
                    Kratos.Expression.LiteralExpressionIO.SetData(exp, 0.0)

                for exp in expressions:
                    for i, value_i in enumerate(self.list_of_discrete_values):
                        if self.residual_type == "exact":
                            partial_gradient_exp = exp.Clone()
                            Kratos.Expression.LiteralExpressionIO.SetData(partial_gradient_exp, 1.0)
                            for j, value_j in enumerate(self.list_of_discrete_values):
                                if i == j:
                                    partial_gradient_exp *= (values - value_j) * 2.0
                                else:
                                    partial_gradient_exp *= (values - value_j) ** 2
                            exp.SetExpression(exp.GetExpression() + partial_gradient_exp.GetExpression())
                        elif self.residual_type == "logarithm":
                            exp.SetExpression(exp.GetExpression() + (((values - value_i) ** (-2)) * (values - value_i) * 2.0).GetExpression())
                    exp.SetExpression(Kratos.Expression.Utils.Collapse(exp).GetExpression())
            else:
                raise RuntimeError(f"Unsupported sensitivity w.r.t. {physical_variable.Name()} requested. Followings are supported sensitivity variables:\n\tSHAPE")
