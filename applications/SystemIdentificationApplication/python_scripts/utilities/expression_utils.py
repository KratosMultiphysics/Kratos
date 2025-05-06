import typing
from enum import Enum
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

class PropertiesDataLocation(Enum):
    ElementProperties = 100
    ConditionProperties = 200

class ExpressionDataLocation(Kratos.Globals.DataLocation):
    ElementProperties = PropertiesDataLocation.ElementProperties,
    ConditionProperties = PropertiesDataLocation.ConditionProperties

ExpressionUnionType = typing.Union[
                            Kratos.Expression.NodalExpression,
                            Kratos.Expression.ConditionExpression,
                            Kratos.Expression.ElementExpression]

def GetContainerExpressionType(data_location: ExpressionDataLocation) -> 'typing.Union[typing.Type[Kratos.Expression.NodalExpression], typing.Type[Kratos.Expression.ConditionExpression], typing.Type[Kratos.Expression.ElementExpression]]':
    if data_location in [ExpressionDataLocation.NodeHistorical, ExpressionDataLocation.NodeNonHistorical]:
        return Kratos.Expression.NodalExpression
    elif data_location == ExpressionDataLocation.Condition or data_location == ExpressionDataLocation.ConditionProperties:
        return Kratos.Expression.ConditionExpression
    elif data_location == ExpressionDataLocation.Element or data_location == ExpressionDataLocation.ElementProperties:
        return Kratos.Expression.ElementExpression
    else:
        raise RuntimeError(f"Unsupported {data_location}.")

def GetContainerExpression(model_part: Kratos.ModelPart, data_location: ExpressionDataLocation, variable: typing.Any) -> ExpressionUnionType:
    cexp_type = GetContainerExpressionType(data_location)
    expression = cexp_type(model_part)
    if data_location == ExpressionDataLocation.NodeHistorical:
        Kratos.Expression.VariableExpressionIO.Read(expression, variable, True)
    elif data_location == ExpressionDataLocation.NodeNonHistorical:
        Kratos.Expression.VariableExpressionIO.Read(expression, variable, False)
    elif data_location == ExpressionDataLocation.Condition:
        Kratos.Expression.VariableExpressionIO.Read(expression, variable)
    elif data_location == ExpressionDataLocation.Element:
        Kratos.Expression.VariableExpressionIO.Read(expression, variable)
    elif data_location == ExpressionDataLocation.ConditionProperties or data_location == ExpressionDataLocation.ElementProperties:
        KratosOA.PropertiesVariableExpressionIO.Read(expression, variable)
    else:
        raise RuntimeError(f"Unsupported {data_location}.")
    return expression

class IntervalBounder:
    def __init__(self, bounds: list[float]) -> None:
        if len(bounds) != 2:
            raise RuntimeError(f"The bounds should be of size 2. [bounds = {bounds}]")
        self.bounds = sorted(bounds)

    def GetBoundGap(self) -> float:
        return self.bounds[1] - self.bounds[0]

    def GetBoundedExpression(self, unbounded_expression: ExpressionUnionType) -> ExpressionUnionType:
        return (unbounded_expression - self.bounds[0]) / self.GetBoundGap()

    def GetUnboundedExpression(self, bounded_expression: ExpressionUnionType) -> ExpressionUnionType:
        return bounded_expression * self.GetBoundGap() + self.bounds[0]
