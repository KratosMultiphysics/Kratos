import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes

def Factory(model_part: Kratos.ModelPart, container_type: Kratos.Globals.DataLocation, filter_radius_settings: Kratos.Parameters) -> ContainerExpressionTypes:
    if not filter_radius_settings.Has("filter_radius_type"):
        raise RuntimeError(f"\"filter_radius_type\" not specified in the following settings:\n{filter_radius_settings}")

    filter_radius_type = filter_radius_settings["filter_radius_type"].GetString()
    if filter_radius_type == "constant":
        defaults = Kratos.Parameters("""{
            "filter_radius_type": "constant",
            "filter_radius"     : 0.2
        }""")
        filter_radius_settings.ValidateAndAssignDefaults(defaults)
        filter_radius = GetContainerExpression(model_part, container_type)
        Kratos.Expression.LiteralExpressionIO.SetData(filter_radius, filter_radius_settings["filter_radius"].GetDouble())
        return filter_radius
    else:
        raise RuntimeError(f"Unsupported filter_radius_type = \"{filter_radius_type}\".")

def GetContainerExpression(model_part: Kratos.ModelPart, container_type: Kratos.Globals.DataLocation) -> ContainerExpressionTypes:
    if container_type == Kratos.Globals.DataLocation.NodeHistorical or container_type == Kratos.Globals.DataLocation.NodeNonHistorical:
        return Kratos.Expression.NodalExpression(model_part)
    elif container_type == Kratos.Globals.DataLocation.Condition:
        return Kratos.Expression.ConditionExpression(model_part)
    elif container_type == Kratos.Globals.DataLocation.Element:
        return Kratos.Expression.ElementExpression(model_part)
    else:
        raise RuntimeError(f"Unsupported container_type.")
