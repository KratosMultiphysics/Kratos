import typing
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.DigitalTwinApplication as KratosDT

def GetNodalSpecificationsWithDirection(model: Kratos.Model, list_of_parameters: 'list[Kratos.Parameters]') -> 'list[KratosDT.Sensors.SensorSpecification]':
    default_settings = Kratos.Parameters("""{
        "name"           : "SENSOR_NAME",
        "model_part_name": "PLEASE_SPECIFY_MODEL_PART_NAME",
        "direction"      : [0.0, 0.0, 0.0],
        "node_ids"       : [],
        "weight"         : 1.0,
        "use_all_nodes"  : false
    }""")

    list_of_specifications: 'list[KratosDT.Sensors.SensorSpecification]' = []
    for parameters in list_of_parameters:
        parameters.ValidateAndAssignDefaults(default_settings)

        sensor_name = parameters["name"].GetString()
        model_part = model[parameters["model_part_name"].GetString()]
        direction = Kratos.Array3(parameters["direction"].GetVector())
        weight = parameters["weight"].GetDouble()

        if parameters["use_all_nodes"].GetBool():
            nodes: 'list[Kratos.Node]' = model_part.Nodes
        else:
            nodes: 'list[Kratos.Node]' = [model_part.GetNode(v) for v in parameters["node_ids"].GetVector()]

        for i, node in enumerate(nodes):
            specification = KratosDT.Sensors.NodalSensorSpecification(sensor_name, i + 1, weight, node)
            specification.SetValue(KratosDT.SENSOR_DIRECTION, direction)
            list_of_specifications.append(specification)

    return list_of_specifications

def PrintSpecificationDataToCSV(list_of_specifications: 'list[KratosDT.Sensors.SensorSpecification]', output_file_name: Path) -> None:
    number_of_clusters = len(list_of_specifications)

    # do nothing if number of clusters is zero
    if number_of_clusters == 0:
        return

    output_file_name.parent.mkdir(exist_ok=True, parents=True)

    with open(str(output_file_name), "w") as csv_output:
        # write the headers of spec
        csv_output.write("#; name; value; location_x; location_y; location_z")
        # now write the headers of data containers
        list_of_vars = []
        for var_name in list_of_specifications[0].GetDataVariableNames():
            csv_output.write(f"; {var_name}")
            list_of_vars.append(Kratos.KratosGlobals.GetVariable(var_name))
        csv_output.write("\n")

        # now write the data
        sensor_id = 0
        for spec in list_of_specifications:
            sensor_id += 1
            loc = spec.GetLocation()
            csv_output.write(f"{sensor_id}; {spec.GetName()}; {spec.GetSensorValue()}; {loc[0]}; {loc[1]}; {loc[2]}")
            for var in list_of_vars:
                csv_output.write(f"; {spec.GetValue(var)}")
            csv_output.write("\n")

def GetNormalizedSpecifications(list_of_specifications: 'list[KratosDT.Sensors.SensorSpecification]', retrieval_mathod: 'typing.Callable[[KratosDT.Sensors.SensorSpecification], typing.Union[Kratos.Expression.NodalExpression, Kratos.Expression.ConditionExpression, Kratos.Expression.ElementExpression]]') -> 'list[KratosDT.Sensors.SensorSpecification]':
    result: 'list[KratosDT.Sensors.SensorSpecification]' = []
    for spec in list_of_specifications:
        cexp = retrieval_mathod(spec)
        norm_l2 = KratosOA.ExpressionUtils.NormL2(cexp)
        if norm_l2 > 0.0:
            cexp.SetExpression(cexp.GetExpression() / norm_l2)
            cexp.SetExpression(cexp.Flatten().GetExpression())
            result.append(spec)
    return result

def GetCosineDistances(list_of_specifications: 'list[KratosDT.Sensors.SensorSpecification]', retrieval_mathod: 'typing.Callable[[KratosDT.Sensors.SensorSpecification], typing.Union[Kratos.Expression.NodalExpression, Kratos.Expression.ConditionExpression, Kratos.Expression.ElementExpression]]') -> 'list[float]':
    """Get the cosine distances in a flat array.

    The distances are computed as follows

         d(u, v) = 1.0 - u\cdot v

    The vectors in u, v must be unit vectors

    Args:
        list_of_specifications (list[KratosDT.Sensors.SensorSpecification]): List of specifications

    Returns:
        list[float]: d(u, v) calculated for all u and v
    """
    results: 'list[float]' = []
    for i, spec_i in enumerate(list_of_specifications):
        cexp_i = retrieval_mathod(spec_i)
        for spec_j in list_of_specifications[i+1:]:
            results.append(1.0 - KratosOA.ExpressionUtils.InnerProduct(cexp_i, retrieval_mathod(spec_j)))
    return results



