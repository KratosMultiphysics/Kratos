import KratosMultiphysics as Kratos
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