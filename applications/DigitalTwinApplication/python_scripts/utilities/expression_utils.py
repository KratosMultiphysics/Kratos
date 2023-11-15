import typing
import KratosMultiphysics as Kratos

ExpressionUnionType = typing.Union[
                            Kratos.Expression.NodalExpression,
                            Kratos.Expression.ConditionExpression,
                            Kratos.Expression.ElementExpression]

def GetContainerExpressionType(data_location: Kratos.Globals.DataLocation) -> 'typing.Union[typing.Type[Kratos.Expression.NodalExpression], typing.Type[Kratos.Expression.ConditionExpression], typing.Type[Kratos.Expression.ElementExpression]]':
    if data_location in [Kratos.Globals.DataLocation.NodeHistorical, Kratos.Globals.DataLocation.NodeNonHistorical]:
        return Kratos.Expression.NodalExpression
    elif data_location == Kratos.Globals.DataLocation.Condition:
        return Kratos.Expression.ConditionExpression
    elif data_location == Kratos.Globals.DataLocation.Element:
        return Kratos.Expression.ElementExpression
    else:
        raise RuntimeError(f"Unsupported {data_location}.")

def GetContainerExpression(model_part: Kratos.ModelPart, data_location: Kratos.Globals.DataLocation, variable: typing.Any) -> ExpressionUnionType:
    cexp_type = GetContainerExpressionType(data_location)
    expression = cexp_type(model_part)
    if data_location == Kratos.Globals.DataLocation.NodeHistorical:
        Kratos.Expression.VariableExpressionIO.Read(expression, variable, True)
    elif data_location == Kratos.Globals.DataLocation.NodeNonHistorical:
        Kratos.Expression.VariableExpressionIO.Read(expression, variable, False)
    elif data_location == Kratos.Globals.DataLocation.Condition:
        Kratos.Expression.VariableExpressionIO.Read(expression, variable)
    elif data_location == Kratos.Globals.DataLocation.Element:
        Kratos.Expression.VariableExpressionIO.Read(expression, variable)
    else:
        raise RuntimeError(f"Unsupported {data_location}.")
    return expression