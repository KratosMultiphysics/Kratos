import typing
import json
from pathlib import Path
from math import sqrt
import numpy as np

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SensorViewUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SupportedValueUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SupportedVariableUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import GetParameterToKratosValuesConverter
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import GetKratosValueToCSVStringConverter
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import GetKratosValueToPythonValueConverter
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import GetNameToCSVString
from KratosMultiphysics.DigitalTwinApplication.utilities.expression_utils import ExpressionUnionType, ExpressionFilterUnionType, GetContainerExpressionType

def GetSensors(model_part: Kratos.ModelPart, list_of_parameters: 'list[Kratos.Parameters]') -> 'list[KratosDT.Sensors.Sensor]':
    """Get list of sensors from given parameters.

    This generates list of sensors from given list of parameters. Each parameter
    corresponds to one settings set for one sensor.

    Args:
        model_part (Kratos.ModelPart): Model part on which the sensors will be generated.
        list_of_parameters (list[Kratos.Parameters]): List of parameters.

    Returns:
        list[KratosDT.Sensors.Sensor]: List of sensors generated.
    """
    point_locator = Kratos.BruteForcePointLocator(model_part)

    defaults_map:'dict[str, Kratos.Parameters]' = {
            "displacement_sensor": Kratos.Parameters(
                                        """{
                                                "type"         : "displacement_sensor",
                                                "name"         : "",
                                                "value"        : 0,
                                                "location"     : [0.0, 0.0, 0.0],
                                                "direction"    : [0.0, 0.0, 0.0],
                                                "weight"       : 1.0,
                                                "variable_data": {}
                                        }""")
    }

    list_of_sensors: 'list[KratosDT.Sensors.Sensor]' = []
    shape_funcs = Kratos.Vector()
    for parameters in list_of_parameters:
        if not parameters.Has("type"):
            raise RuntimeError(f"The sensor parameters does not contain \"type\".")
        sensor_type = parameters["type"].GetString()

        if not sensor_type in defaults_map.keys():
            raise RuntimeError(f"Unsupported sensor type = \"{sensor_type}\" requested. Followings are supported:\n\t" + "\n\t".join(defaults_map.keys()))
        parameters.ValidateAndAssignDefaults(defaults_map[sensor_type])

        if sensor_type == "displacement_sensor":
            name = parameters["name"].GetString()
            loc = parameters["location"].GetVector()
            loc = Kratos.Point(loc[0], loc[1], loc[2])
            weight = parameters["weight"].GetDouble()
            direction = parameters["direction"].GetVector()
            elem_id = point_locator.FindElement(loc, shape_funcs, Kratos.Configuration.Initial, 1e-8)
            sensor = KratosDT.Sensors.DisplacementSensor(name, loc, direction, model_part.GetElement(elem_id), weight)
            AddSensorVariableData(sensor, parameters["variable_data"])
            list_of_sensors.append(sensor)
    return list_of_sensors

def PrintSensorListToCSV(output_file_name: Path, list_of_sensors: 'list[KratosDT.Sensors.Sensor]', list_of_sensor_properties: 'list[str]') -> None:
    """Writes data in the list_of_sensors to CSV file.

    This method writes data from all the list_of_sensors to the CSV file. The columns of the
    CSV is determined by list_of_sensor_properties.

    Args:
        output_file_name (Path): CSV file name (including the .csv extension).
        list_of_sensors (list[Sensor]): List of sensors.
        list_of_sensor_properties (list[str]): List of columns to write.
    """
    number_of_sensor_views = len(list_of_sensors)

    # do nothing if number of clusters is zero
    if number_of_sensor_views == 0:
        return

    # create output path
    output_file_name.parent.mkdir(exist_ok=True, parents=True)

    with open(str(output_file_name), "w") as file_output:
        # print the headers
        first_sensor = list_of_sensors[0]

        # write the csv id
        file_output.write("#")

        # store converters to CSV in dicts
        params_value_dict: 'dict[str, typing.Callable[[Kratos.Parameters], str]]' = {}
        params_str_dict: 'dict[str, typing.Callable[[SupportedValueUnionType], str]]' = {}
        data_value_dict: 'dict[SupportedVariableUnionType, typing.Callable[[SupportedValueUnionType], str]]' = {}

        dummy_params = first_sensor.GetSensorParameters()
        list_of_params_keys: 'list[str]' = []
        list_of_data_keys: 'list[SupportedVariableUnionType]' = []
        param_property_headers = ""
        value_property_headers = ""
        for sensor_property in list_of_sensor_properties:
            if dummy_params.Has(sensor_property):
                list_of_params_keys.append(sensor_property)
                value_func = GetParameterToKratosValuesConverter(dummy_params[sensor_property])
                value = value_func(dummy_params[sensor_property])
                params_value_dict[sensor_property] = value_func
                params_str_dict[sensor_property] = GetKratosValueToCSVStringConverter(value)
                param_property_headers += f",{GetNameToCSVString(sensor_property, value)}"
            else:
                var = Kratos.KratosGlobals.GetVariable(sensor_property)
                if first_sensor.Has(var):
                    list_of_data_keys.append(var)
                    value = first_sensor.GetValue(var)
                    data_value_dict[var] = GetKratosValueToCSVStringConverter(value)
                    value_property_headers += f",{GetNameToCSVString(sensor_property, value)}"
                else:
                    raise RuntimeError(f"Sensor property = \"{sensor_property}\" not found in either sensor parameters nor in data value container. [ Sensor = {first_sensor} ].")

        file_output.write(f"{param_property_headers}{value_property_headers}\n")

        for i, sensor in enumerate(list_of_sensors):
            # write the id
            file_output.write(f"{i+1}")

            # first write the values from sensor params
            sensor_params = sensor.GetSensorParameters()
            for sensor_property in list_of_params_keys:
                file_output.write(f",{params_str_dict[sensor_property](params_value_dict[sensor_property](sensor_params[sensor_property]))}")

            # now write the data values from sensor
            for var in list_of_data_keys:
                file_output.write(f",{data_value_dict[var](sensor.GetValue(var))}")

            file_output.write("\n")

def GetNormalizedSensorViews(list_of_sensor_views: 'list[SensorViewUnionType]') -> 'list[SensorViewUnionType]':
    """Get the normalized sensor views.

    This method returns sensor views which will contain normalized container expressions.
    The sensor views with zero norm will not be included in the return list

    Args:
        list_of_sensor_views (list[SensorViewUnionType]): List of sensors to normalize.

    Returns:
        list[SensorViewUnionType]: Normalized list of sensor views.
    """
    result: 'list[SensorViewUnionType]' = []
    for sensor_view in list_of_sensor_views:
        cexp = sensor_view.GetContainerExpression()
        norm_l2 = KratosOA.ExpressionUtils.NormL2(cexp)
        if norm_l2 > 0.0:
            cexp.SetExpression(cexp.GetExpression() / norm_l2)
            cexp.SetExpression(cexp.Flatten().GetExpression())
            result.append(sensor_view)
    return result

def GetCosineDistances(list_of_sensor_views: 'list[SensorViewUnionType]') -> 'list[float]':
    """Get the cosine distances in a flat array.

    The distances are computed as follows

         d(u, v) = 1.0 - u\cdot v

    The vectors in u, v must be unit vectors.

    Args:
        list_of_sensor_views (list[SensorViewUnionType]): List of sensor views

    Returns:
        list[float]: d(u, v) calculated for all u and v
    """
    results: 'list[float]' = []
    for i, sensor_view_i in enumerate(list_of_sensor_views):
        for sensor_view_j in list_of_sensor_views[i+1:]:
            results.append(1.0 - KratosOA.ExpressionUtils.InnerProduct(sensor_view_i.GetContainerExpression(), sensor_view_j.GetContainerExpression()))
    return results

def GetEuclideanDistances(list_of_sensor_views: 'list[SensorViewUnionType]') -> 'list[float]':
    """Get the Euclidean distances in a flat array.

    The distances are computed as follows

         d(u, v) = |u-v|

    Args:
        list_of_sensor_views (list[SensorViewUnionType]): List of sensor views

    Returns:
        list[float]: d(u, v) calculated for all u and v
    """
    results: 'list[float]' = []
    for i, sensor_view_i in enumerate(list_of_sensor_views):
        for sensor_view_j in list_of_sensor_views[i+1:]:
            dist = sensor_view_i.GetSensor().GetLocation() - sensor_view_j.GetSensor().GetLocation()
            results.append((dist[0]**2 + dist[1]**2 + dist[2]**2) ** 0.5)
    return results

def GetSubDistances(original_distances: 'list[float]', sub_indices_list: 'list[int]') -> 'list[float]':
    m = len(original_distances)
    n = int((1 + sqrt(8*m + 1)) / 2)
    if m != (n * (n-1)) // 2:
        raise RuntimeError("Invalid original distances vector given.")

    sub_indices_list = sorted(sub_indices_list)
    sub_distances: 'list[float]' = []
    for i, sub_index_i in enumerate(sub_indices_list):
        for sub_index_j in sub_indices_list[i+1:]:
            sub_distances.append(original_distances[n * sub_index_i + sub_index_j - ((sub_index_i + 2) * (sub_index_i + 1)) // 2])
    return sub_distances

def AddSensorVariableData(sensor: KratosDT.Sensors.Sensor, variable_data: Kratos.Parameters) -> None:
    """Adds sensor variable data.

    Args:
        sensor (KratosDT.Sensors.Sensor): Sensor to add data.
        variable_data (Kratos.Parameters): Parameters containing variable data.
    """
    for var_name, var_value in variable_data.items():
        var = Kratos.KratosGlobals.GetVariable(var_name)
        value_func =  GetParameterToKratosValuesConverter(var_value)
        sensor.SetValue(var, value_func(var_value))

def GetSensorTypeDict(list_of_sensors: 'list[SensorViewUnionType]') -> 'dict[str, list[KratosDT.Sensors.Sensor]]':
    """Get the sensor type dict from list of sensors with keys being type of the sensor.

    Args:
        list_of_sensors (list[KratosDT.Sensors.Sensor]): List of sensors

    Returns:
        dict[str, list[KratosDT.Sensors.Sensor]]: Sensor type dict.
    """
    result: 'dict[str, list[KratosDT.Sensors.Sensor]]' = {}
    for sensor in list_of_sensors:
        sensor_type = type(sensor.GetSensor()).__name__
        if sensor_type not in result.keys():
            result[sensor_type] = []
        result[sensor_type].append(sensor)
    return result

def PrintSensorListToJson(output_file_name: Path, list_of_sensors: 'list[KratosDT.Sensors.Sensor]') -> None:
    number_of_sensor_views = len(list_of_sensors)

    # do nothing if number of clusters is zero
    if number_of_sensor_views == 0:
        return

    # create output path
    output_file_name.parent.mkdir(exist_ok=True, parents=True)

    list_of_vars: 'list[tuple[SupportedVariableUnionType, typing.Callable[[SupportedValueUnionType], typing.Union[bool, int, float, str, list[float]]]]]' = []
    for var_name in list_of_sensors[0].GetDataVariableNames():
        var: SupportedVariableUnionType = Kratos.KratosGlobals.GetVariable(var_name)
        list_of_vars.append((var, GetKratosValueToPythonValueConverter(list_of_sensors[0].GetValue(var))))

    json_sensors = {"list_of_sensors": []}
    with open(str(output_file_name), "w") as file_output:
        for sensor in list_of_sensors:
            json_params = json.loads(sensor.GetSensorParameters().WriteJsonString())
            json_params["variable_data"] = {}

            for var, func in list_of_vars:
                json_params["variable_data"][var.Name()] = func(sensor.GetValue(var))
            json_sensors["list_of_sensors"].append(json_params)

        file_output.write(json.dumps(json_sensors, indent=4))

def GetBestCoverageSensorView(list_of_sensor_views: 'list[SensorViewUnionType]') -> SensorViewUnionType:
    if len(list_of_sensor_views) == 0:
        raise RuntimeError("No sensor views found.")

    front_sensor_view = list_of_sensor_views[0]
    overall_updating_exp = front_sensor_view.GetContainerExpression().Clone()
    overall_updating_exp = overall_updating_exp * 0.0 + 1.0
    overall_updating_exp /= KratosOA.ExpressionUtils.NormL2(overall_updating_exp)
    overall_updating_exp.SetExpression(overall_updating_exp.Flatten().GetExpression())
    return min(list_of_sensor_views, key=lambda x: 1.0 - KratosOA.ExpressionUtils.InnerProduct(overall_updating_exp, x.GetContainerExpression()))

def GetDistance(i: int, j: int, compressed_distance_matrix: 'list[float]') -> float:
    local_i_j = sorted([i, j])
    n = len(compressed_distance_matrix)
    m = int((1.0 + sqrt(8.0*n + 1.0)) / 2)
    return compressed_distance_matrix[m*local_i_j[0] + local_i_j[1] - ((local_i_j[0] + 2)*(local_i_j[0] + 1)) // 2]

def ComputeHeatMap(list_of_sensor_views: 'list[SensorViewUnionType]') -> ExpressionUnionType:
    if len(list_of_sensor_views) == 0:
        raise RuntimeError("No sensor views found.")

    front_sensor_view = list_of_sensor_views[0]
    heat_map = front_sensor_view.GetContainerExpression().Clone()
    for sensor_view in list_of_sensor_views[1:]:
        heat_map += sensor_view.GetContainerExpression()
    heat_map /= KratosOA.ExpressionUtils.NormL2(heat_map)
    return heat_map.Flatten()

def GetFilter(model_part: Kratos.ModelPart, filter_field_location: Kratos.Globals.DataLocation, filter_parameters: Kratos.Parameters) -> ExpressionFilterUnionType:
    if not filter_parameters.Has("filter_type"):
        raise RuntimeError(f"Parameters does not contain filter_type. {filter_parameters}")

    filter_type = filter_parameters["filter_type"].GetString()
    if filter_type == "entity_nodal_entity_filter":
        if filter_field_location == Kratos.Globals.DataLocation.Condition:
            return KratosOA.ConditionNodeConditionFilter(model_part)
        elif filter_field_location == Kratos.Globals.DataLocation.Element:
            return KratosOA.ElementNodeElementFilter(model_part)
        else:
            raise RuntimeError(f"Unsupported filter field location = {filter_field_location.name}")
    elif filter_type == "explicit_vertex_morphing":
        default_settings = Kratos.Parameters("""{
            "filter_type"               : "explicit_vertex_morphing",
            "filter_radius"             : 5.0,
            "filter_function_type"      : "linear",
            "fixed_model_part_name"     : "",
            "damping_function_type"     : "sigmoidal",
            "max_nodes_in_filter_radius": 1000
        }""")
        filter_parameters.ValidateAndAssignDefaults(default_settings)

        expression_type = GetContainerExpressionType(filter_field_location)
        filter_radius_exp = expression_type(model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(filter_radius_exp, filter_parameters["filter_radius"].GetDouble())
        fixed_model_part_name = filter_parameters["fixed_model_part_name"].GetString()
        filter_function_type = filter_parameters["filter_function_type"].GetString()
        damping_function_type = filter_parameters["damping_function_type"].GetString()
        max_filtering_nodes = filter_parameters["max_nodes_in_filter_radius"].GetInt()

        if expression_type == Kratos.Expression.NodalExpression:
            filter_type = KratosOA.NodalExplicitFilter
        elif expression_type == Kratos.Expression.ConditionExpression:
            filter_type = KratosOA.ConditionExplicitFilter
        elif expression_type == Kratos.Expression.ElementExpression:
            filter_type = KratosOA.ElementExplicitFilter
        else:
            raise RuntimeError("Unsupported filter type.")

        if fixed_model_part_name == "":
            vm_filter = filter_type(model_part, filter_function_type, max_filtering_nodes)
        else:
            vm_filter = filter_type(model_part, model_part.GetModel()[fixed_model_part_name], filter_function_type, damping_function_type, max_filtering_nodes)

        vm_filter.SetFilterRadius(filter_radius_exp)
        return vm_filter
    else:
        raise RuntimeError(f"Unsupported filter_type = {filter_type}.")

def GetSimilarSensorViews(reference_expression: ExpressionUnionType, cosine_similarity_lb: float, cosine_similarity_ub: float, list_of_sensor_views: 'list[SensorViewUnionType]') -> 'list[SensorViewUnionType]':
    results: 'list[SensorViewUnionType]' = []
    for sensor_view in list_of_sensor_views:
        cosine_similarity = KratosOA.ExpressionUtils.InnerProduct(sensor_view.GetContainerExpression(), reference_expression)
        if cosine_similarity <= cosine_similarity_ub and cosine_similarity > cosine_similarity_lb:
            results.append(sensor_view)
    return results

def GetSensorCoverageMasks(sensor_views: 'list[SensorViewUnionType]', domain_size_exp: ExpressionUnionType) -> 'list[ExpressionUnionType]':
    list_of_coverage_masks: 'list[ExpressionUnionType]' = []
    total_domain_size = KratosDT.SensorUtils.Sum(domain_size_exp)

    for sensor_view in sensor_views:
        _, coverage_mask = KratosDT.SensorUtils.GetEntityCoverageMask(sensor_view)
        coverage = KratosDT.SensorUtils.Sum(domain_size_exp.Scale(coverage_mask)) / total_domain_size
        Kratos.Logger.PrintInfo("GetSensorCoverageMasks", f"Coverage of \"{sensor_view.GetSensor().GetName()}\" is {coverage * 100:6.3f}%")
        sensor_view.GetSensor().SetValue(KratosDT.SENSOR_COVERAGE, coverage)
        list_of_coverage_masks.append(coverage_mask)
    return list_of_coverage_masks

def ComputeMinimumDistanceSquare(reference_view: SensorViewUnionType, sensor_views: 'list[SensorViewUnionType]') -> float:
    minimum_distance = 1e+100
    for sensor_view in sensor_views:
        current_distance = reference_view.GetSensor().GetLocation() - sensor_view.GetSensor().GetLocation()
        current_distance = current_distance[0] ** 2 + current_distance[1] ** 2 + current_distance[2] ** 2
        if minimum_distance > current_distance:
            minimum_distance = current_distance
    return minimum_distance

def GetMostDistancedMax(relaxation: float, list_of_values: 'list[float]', potential_sensor_views: 'list[SensorViewUnionType]', current_sensor_views: 'list[SensorViewUnionType]') -> int:
    list_of_values_with_index = [(i, v) for i, v in enumerate(list_of_values)]
    list_of_values_with_index = sorted(list_of_values_with_index, key=lambda x: x[1], reverse=True)

    cut_off_value = list_of_values_with_index[0][1] * (1 - relaxation)

    max_distance = 0.0
    max_distanced_index = 0
    if len(current_sensor_views) > 0:
        for i, v in list_of_values_with_index:
            if v >= cut_off_value:
                min_distance = 1e+9
                for sensor_view in current_sensor_views:
                    current_distance = sensor_view.GetSensor().GetLocation() - potential_sensor_views[i].GetSensor().GetLocation()
                    current_distance = current_distance[0] ** 2 + current_distance[1] ** 2 + current_distance[2] ** 2
                    min_distance = min(min_distance, current_distance)
                if max_distance < min_distance:
                    max_distance = min_distance
                    max_distanced_index = i
            else:
                break
    return max_distanced_index

def GetMostDistancedMin(relaxation: float, list_of_values: 'list[float]', potential_sensor_views: 'list[SensorViewUnionType]', current_sensor_views: 'list[SensorViewUnionType]') -> int:
    list_of_values_with_index = [(i, v) for i, v in enumerate(list_of_values)]
    list_of_values_with_index = sorted(list_of_values_with_index, key=lambda x: x[1])

    cut_off_value = list_of_values_with_index[0][1] * (1 + relaxation)

    max_distance = 0.0
    max_distanced_index = 0
    if len(current_sensor_views) > 0:
        for i, v in list_of_values_with_index:
            if v <= cut_off_value:
                min_distance = 1e+9
                for sensor_view in current_sensor_views:
                    current_distance = sensor_view.GetSensor().GetLocation() - potential_sensor_views[i].GetSensor().GetLocation()
                    current_distance = current_distance[0] ** 2 + current_distance[1] ** 2 + current_distance[2] ** 2
                    min_distance = min(min_distance, current_distance)
                if max_distance < min_distance:
                    max_distance = min_distance
                    max_distanced_index = i
            else:
                break
    return max_distanced_index

